import time
import json
import ctypes
import subprocess
import numpy as np
import comtypes.client

from pathlib import Path
from typing import Protocol
from datetime import datetime
from comtypes import automation
from comtypes.automation import VT_R8
from pythoncom import CoInitialize, CoUninitialize


class OpenSTAADProperty(Protocol):
    pass
class OpenSTAADGeometry(Protocol):
    pass
class OpenSTAADSupport(Protocol):
    pass
class OpenSTAADLoad(Protocol):
    pass
class OpenSTAADCommand(Protocol):
    pass
class OpenSTAADOutput(Protocol):
    pass


class OpenSTAAD(Protocol):
    @property
    def Property(self) -> OpenSTAADProperty: ...
    @property
    def Geometry(self) -> OpenSTAADGeometry: ...
    @property
    def Support(self) -> OpenSTAADSupport: ...
    @property
    def Load(self) -> OpenSTAADLoad: ...
    @property
    def Command(self) -> OpenSTAADCommand: ...
    @property
    def Output(self) -> OpenSTAADOutput: ...


def make_variant_vt_ref(obj: any, var_type: int) -> automation.VARIANT:
    """Wraps an object in a VARIANT with VT_BYREF flag."""
    var = automation.VARIANT()
    var._.c_void_p = ctypes.addressof(obj)
    var.vt = var_type | automation.VT_BYREF
    return var


def make_safe_array_long(size: int):
    """Creates a SAFEARRAY of longs of given size."""
    return automation._midlSAFEARRAY(ctypes.c_long).create([0] * size)


def make_safe_array_double(size: int):
    """Creates a SAFEARRAY of doubles of given size."""
    return automation._midlSAFEARRAY(ctypes.c_double).create([0] * size)


def make_safe_array_double_values(values: list[float]):
    """Creates a SAFEARRAY of doubles from a list of values."""
    return automation._midlSAFEARRAY(ctypes.c_double).create(values)


def make_safe_array_long_values(values: list[int]):
    """Create a SAFEARRAY of longs from a list of values."""
    return automation._midlSAFEARRAY(ctypes.c_long).create(values)


def assign_shell_thickness(staad_property: OpenSTAADProperty, thickness_list: list[float]) -> int:
    """Assigns plate thickness property via STAAD API."""
    thickness_values = np.full(len(thickness_list), thickness_list[0], dtype=np.double)
    thickness_list = thickness_values.tolist()
    safe_array_thickness = make_safe_array_double_values(thickness_list)
    variant_thickness = make_variant_vt_ref(safe_array_thickness, automation.VT_ARRAY | automation.VT_R8)
    staad_property._FlagAsMethod("CreatePlateThicknessProperty")
    return staad_property.CreatePlateThicknessProperty(variant_thickness)


def rigid_mp_constraints(
    staad_property: OpenSTAADProperty, control_node: int, dependent_nodes: list[int] | np.ndarray) -> int:
    """Calls STAAD API's AddControlDependentRelation method."""
    if isinstance(dependent_nodes, np.ndarray):
        dependent_nodes = dependent_nodes.tolist()
    safe_array_dependent_nodes = make_safe_array_long_values(dependent_nodes)
    variant_dependent_nodes = make_variant_vt_ref(safe_array_dependent_nodes, automation.VT_ARRAY | automation.VT_I4)
    staad_property._FlagAsMethod("AddControlDependentRelation")
    return staad_property.AddControlDependentRelation(
        control_node, -1, 1, 1, 1, 1, 1, 1, variant_dependent_nodes
    )


def get_all_hole_node_tags(hole_boundaries: dict[str, list[int]]) -> list[int]:
    """Returns all hole node tags (each incremented by 1 due to the geometry
    generation functions in models)."""
    return [hole + 1 for holes in hole_boundaries.values() for hole in holes]


def convert_vkt_dir_2_staad(l_dir: str) -> str:
    """Converts Viktor direction to STAAD equivalent."""
    conversion_dict = {
        "Y+": "Z+",
        "Y-": "Z-",
        "Z+": "Y+",
        "Z-": "Y-"
    }
    return conversion_dict.get(l_dir, l_dir)


def create_point_loads(load,case_num: int, load_list: list, top_cylinder_nodetag: int)->int:
    ret = load.SetLoadActive(case_num)
    for args in load_list:
        ld_dir = convert_vkt_dir_2_staad(args["ld_dir"])
        mag = args["mag"]
        # Point Loads:
        if args["type"] == "Point Load":
            if ld_dir == "Y+":
                load.AddNodalLoad(top_cylinder_nodetag, 0, mag, 0, 0, 0, 0)
            elif ld_dir == "Y-":
                load.AddNodalLoad(top_cylinder_nodetag, 0, -mag, 0, 0, 0, 0)
            elif ld_dir == "Z+":
                load.AddNodalLoad(top_cylinder_nodetag, 0, 0, mag, 0, 0, 0)
            elif ld_dir == "Z-":
                load.AddNodalLoad(top_cylinder_nodetag, 0, 0, -mag, 0, 0, 0)
            elif ld_dir == "X+":
                load.AddNodalLoad(top_cylinder_nodetag, mag, 0, 0, 0, 0, 0)
            elif ld_dir == "X-":
                load.AddNodalLoad(top_cylinder_nodetag, -mag, 0, 0, 0, 0, 0)
        else:
            # Moment Load in kN*mm then mag*1000 converts it to kN*m
            if ld_dir == "Y+":
                load.AddNodalLoad(top_cylinder_nodetag, 0, 0, 0, 0, mag*1000, 0)
            elif ld_dir == "Y-":
                load.AddNodalLoad(top_cylinder_nodetag, 0, 0, 0, 0, -mag*1000, 0)
            elif ld_dir == "Z+":
                load.AddNodalLoad(top_cylinder_nodetag, 0, 0, 0, 0, 0, mag*1000)
            elif ld_dir == "Z-":
                load.AddNodalLoad(top_cylinder_nodetag, 0, 0, 0, 0, 0, -mag*1000)
            elif ld_dir == "X+":
                load.AddNodalLoad(top_cylinder_nodetag, 0, 0, 0, mag*1000, 0, 0)
            elif ld_dir == "X-":
                load.AddNodalLoad(top_cylinder_nodetag, 0, 0, 0, -mag*1000, 0, 0)
    return ret


def base_plate_compression_only_supports(
                                  triangles: dict[str, dict], 
                                  nodes_w_constraints: list[int], 
                                  hole_boundaries: dict[str, list[int]], 
                                  staad_support: OpenSTAADSupport, 
                                  support_id: int):

    """Assigns support to base plate entities via STAAD API."""
    # Collect entity IDs from triangles where the group is "base plate"
    entity_list = [int(tid) for tid, tdata in triangles.items() if tdata.get("group") == "base plate"]
    # Wrap the support_id in a VARIANT.
    variant_support_id = automation.VARIANT(support_id)
    variant_support_id.vt = automation.VT_I4
    # If the entity_list happens to be a NumPy array, convert it to a list.
    if isinstance(entity_list, np.ndarray):
        entity_list = entity_list.tolist()
    # Create a SAFEARRAY for the entity list using c_long.
    safe_array_entity = make_safe_array_long_values(entity_list)
    # Wrap the safe array in a VARIANT.
    variant_entity = make_variant_vt_ref(safe_array_entity, automation.VT_ARRAY | automation.VT_I4)
    # Flag the method if required by the STAAD API.
    staad_support._FlagAsMethod("AssignSupportToEntityList")
    # Call the method with the VARIANT parameters.
    return staad_support.AssignSupportToEntityList(variant_support_id, variant_entity)


def get_displacements(output: OpenSTAADOutput, node_i: int, lc: int = 1) -> dict[str, dict[str, float]]:
    """Retrieves node displacements (converted to mm)."""
    safe_n1 = make_safe_array_double(6)
    x = make_variant_vt_ref(safe_n1, automation.VT_ARRAY | automation.VT_R8)
    output.GetNodeDisplacements(node_i, lc, x)
    return {str(node_i): {"x": x.value[0][0] * 25.4, "y": x.value[0][1] * 25.4, "z": x.value[0][2] * 25.4}}


def get_reaction_loads(output: OpenSTAADOutput, node_id: int, lc: int = 1) -> dict[str, dict[str, float]]:
    """Retrieves support reactions for a node."""
    safe_n1 = make_safe_array_double(6)
    x = make_variant_vt_ref(safe_n1,  automation.VT_ARRAY |  automation.VT_R8)
    output.GetSupportReactions(int(node_id),lc,x)
    Fx, Fy, Fz, *_  = x.value[0]
    return {str(node_id): {"Fx [kN]": round(Fx, 4), "Fy [kN]": round(Fy, 2), "Fz [kN]": round(Fz, 2)}}


def get_all_reactions(output: OpenSTAADOutput, nodes_w_constraints: list[int], lc: int = 1) -> dict[str, dict[str, float]]:
    reactions_dict: dict[str,dict[str,float]] = {}
    for bolt in nodes_w_constraints[:-1]:
        reactions_dict.update(get_reaction_loads(output=output, node_id=bolt, lc=lc))
    return reactions_dict


def get_all_displacements(output, triangles: dict, lc: int = 1):
    disp_dict: dict[str, dict] = {}
    for _, triangles_args in triangles.items():
        node_i_tag = triangles_args["NodeI"]
        node_j_tag = triangles_args["NodeJ"]
        node_k_tag = triangles_args["NodeK"]

        for node_tag in [node_i_tag, node_j_tag, node_k_tag]:
            if str(node_tag) not in disp_dict:
                disp_dict.update(get_displacements(output=output, node_i=int(node_tag), lc=lc))
    return disp_dict


def rigid_links(staad_property: OpenSTAADProperty)->int:
    staad_property._FlagAsMethod("CreateIsotropicMaterialProperties")
    staad_property._FlagAsMethod("CreatePrismaticCircleProperty")
    material_name = "RigidMaterial2"
    staad_property.CreateIsotropicMaterialProperties(
        material_name,
        1e15,  # E
        0.98,  # u
        0.06,  # Gamma
        2.78e-06,  # Thermal
        0.06,  # Damping
    )
    staad_property.SetMaterialName(material_name)
    # Create a prismatic circle property with huge diameter to create a "rigid" element
    return staad_property.CreatePrismaticCircleProperty(10000000)


def create_rigid_links(
    staad_property: OpenSTAADProperty, geometry: OpenSTAADGeometry, nodes: dict[str, dict], top_nodes_cylinder: dict[str, list[int]], nodes_w_constraints: list[int], prop_num: int
) -> None:
    """Creates rigid beam links between nodes."""
    rigid_links = []
    master_node = nodes_w_constraints[-1]
    cords = nodes[str(master_node)]
    geometry.CreateNode(int(master_node), cords["x"], cords["z"], cords["y"])
    for node_tag in top_nodes_cylinder["top"][:-1]:
        beam_tag = geometry.AddBeam(int(node_tag) + 1, int(master_node))
        rigid_links.append(beam_tag)
    [staad_property.AssignBeamProperty(int(line_id), prop_num) for line_id in rigid_links]


def create_top_contraints(staad_property: OpenSTAADProperty, hole_boundaries: dict[str,list[int]], master: int)->None:
    """ Creates Master Slave relationships between the nodes - This dont work if you the master node
    or depended nodes have a constraint like a support"""
    dependent_nodes = [node_tag+1 for node_tag in hole_boundaries["top"][:-1]]
    new_master = dependent_nodes[-1]
    dependent_nodes[-1] = master
    rigid_mp_constraints(staad_property,new_master, dependent_nodes)


def create_bolt_contraints(
    staad_property: OpenSTAADProperty, geometry: OpenSTAADGeometry, nodes: dict[str, dict], hole_boundaries: dict[str, list[int]], nodes_w_constraints: list[int], prop_num: int
) -> None:
    """Creates bolt constraint with rigid links between nodes."""
    rigid_links = []
    for master_node, hole_node_list in zip(nodes_w_constraints[:-1], hole_boundaries.values()):
        cords = nodes[str(master_node)]
        geometry.CreateNode(int(master_node), cords["x"], cords["z"], cords["y"])
        for node_tag in hole_node_list:
            beam_tag = geometry.AddBeam(int(node_tag + 1), int(master_node))
            rigid_links.append(beam_tag)

    [staad_property.AssignBeamProperty(int(line_id), prop_num) for line_id in rigid_links]


def create_bolt_support(support:OpenSTAADSupport)->int:
    """Creates a fixed (but with supports) support using STAAD API."""
    # Define the specifications as lists of 6 double values.
    varReleaseSpec = [0, -1, 0, 0, -1, 0]       # FX, FY, FZ, MX, MY, MZ
    varSpringSpec  = [0, 10000, 0, 0, 1000000000, 0]  # KFX, KFY, KFZ, KMX, KMY, KMZ
    # Create safe arrays for both arrays.
    safe_array_release = make_safe_array_double_values(varReleaseSpec)
    safe_array_spring  = make_safe_array_double_values(varSpringSpec)
    # Wrap the safe arrays in VARIANTs.
    # VT_R8 corresponds to double values.
    variant_release = make_variant_vt_ref(safe_array_release, automation.VT_ARRAY | automation.VT_R8)
    variant_spring  = make_variant_vt_ref(safe_array_spring,  automation.VT_ARRAY | automation.VT_R8)
    support._FlagAsMethod("CreateSupportFixedBut")
    # Call the method with the VARIANT parameters.
    support_id = support.CreateSupportFixedBut(variant_release, variant_spring)
    return support_id


def get_max_von_misses(output: OpenSTAADOutput, plate_id: int, lc: int)->float:
    """Returns maximum von Mises stress for a plate (from kilo pound/square inch to MPa)."""
    # Create cytype
    pdvONT = ctypes.c_double()
    pdvONB = ctypes.c_double()
    # Create a variant reft
    pdVONT = make_variant_vt_ref(pdvONT, VT_R8)
    pdVONB = make_variant_vt_ref(pdvONB, VT_R8)
    # Get the top and bottom Von Misess Stresses
    output.GetPlateCenterVonMisesStresses(plate_id, lc, pdVONT, pdVONB)
    return max([pdvONT.value, pdvONB.value]) * 6.89476 # Convert to Mpa


def get_all_von_misses(output: OpenSTAADOutput, triangles: dict[str, dict], lc: int = 1) -> dict[str, int]:
    max_stress_restuls = {}
    for triangle_id, vals in triangles.items():
        max_stress_restuls[str(triangle_id)] = get_max_von_misses(output=output, plate_id=int(triangle_id), lc=lc)
    return max_stress_restuls


def run_staad()->None:
    CoInitialize()
    # Replace with your version and file path.
    staad_path = r"C:\Program Files\Bentley\Engineering\STAAD.Pro 2024\STAAD\Bentley.Staad.exe"
    # Launch STAAD.Pro
    staad_process = subprocess.Popen([staad_path])
    print("Launching STAAD.Pro...")
    time.sleep(18)
    # Connect to OpenSTAAD.
    openstaad = comtypes.client.GetActiveObject("StaadPro.OpenSTAAD")

    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M")
    std_file_path = Path.cwd() / f"Structure_{timestamp}.std"
    length_unit = 5  # mm.
    force_unit = 5  # Kilo Newton.
    openstaad.NewSTAADFile(str(std_file_path), length_unit, force_unit)

    # Load lines and nodes
    # # Create joints
    input_json = Path.cwd() / "inputs.json"
    with open(input_json) as jsonfile:
        data = json.load(jsonfile)

    nodes, triangles, nodes_w_constraints, hole_boundaries, top_nodes, load_list, plate_thickness = data[:]

    # Wait to load interface
    time.sleep(10)

    # Set Material and Beam Section
    staad_property = openstaad.Property
    staad_property._FlagAsMethod("SetMaterialName")
    staad_property._FlagAsMethod("AssignMaterialToPlate")
    staad_property._FlagAsMethod("AssignBeamProperty")
    material_name = "STEEL"
    staad_property.SetMaterialName(material_name)

    # Create Members.
    geometry = openstaad.Geometry
    geometry._FlagAsMethod("CreateNode")
    geometry._FlagAsMethod("CreatePlate")
    geometry._FlagAsMethod("AddBeam")

    staad_property._FlagAsMethod("AssignPlateThickness")
    staad_property._FlagAsMethod("AssignBeamProperty")
    staad_property._FlagAsMethod("CreatePlateThicknessProperty")
    spec_no = assign_shell_thickness(staad_property=staad_property,thickness_list=[plate_thickness, plate_thickness, plate_thickness, plate_thickness])

    create_nodes = set()
    for triangle_id, vals in triangles.items():
        node_i_id = str(vals["NodeI"])
        node_i_cords = nodes[node_i_id]

        node_j_id = str(vals["NodeJ"])
        node_j_cords = nodes[node_j_id]

        node_k_id = str(vals["NodeK"])
        node_k_cords = nodes[node_k_id]

        if node_i_id not in create_nodes:
            geometry.CreateNode(int(node_i_id), node_i_cords["x"], node_i_cords["z"], node_i_cords["y"])
        if node_j_id not in create_nodes:
            geometry.CreateNode(int(node_j_id), node_j_cords["x"], node_j_cords["z"], node_j_cords["y"])

        if node_k_id not in create_nodes:
            geometry.CreateNode(int(node_k_id), node_k_cords["x"], node_k_cords["z"], node_k_cords["y"])
        ret = geometry.CreatePlate(int(triangle_id), int(node_i_id), int(node_j_id), int(node_k_id))

        # Assign beam property to beam ids.
        ret = staad_property.AssignPlateThickness(int(triangle_id), spec_no)
        ret = staad_property.AssignMaterialToPlate(material_name, int(triangle_id))

    # Create supports.
    staad_support = openstaad.Support
    staad_support._FlagAsMethod("CreatePlateMat")
    staad_support._FlagAsMethod("AssignSupportToNode")
    staad_support._FlagAsMethod("CreateSupportFixedBut")

    elastic_support_id = staad_support.CreatePlateMat(1, 0.03, 0, 1)

    base_plate_compression_only_supports(
        triangles=triangles, nodes_w_constraints=nodes_w_constraints, hole_boundaries=hole_boundaries, staad_support=staad_support, support_id=elastic_support_id
    )

    staad_support._FlagAsMethod("RemoveSupportFromNode")
    hole_nodes = get_all_hole_node_tags(hole_boundaries=hole_boundaries)
    [staad_support.RemoveSupportFromNode(tag) for tag in hole_nodes]

    create_bolt_support(support=staad_support)
    [staad_support.RemoveSupportFromNode(int(node)) for node in nodes_w_constraints[:-1]]

    for node in nodes_w_constraints[:-1]:
        staad_support.AssignSupportToNode(int(node), elastic_support_id + 1)

    prop_num = rigid_links(staad_property=staad_property)
    create_bolt_contraints(
        staad_property=staad_property, geometry=geometry, nodes=nodes, hole_boundaries=hole_boundaries, nodes_w_constraints=nodes_w_constraints, prop_num=prop_num
    )
    create_rigid_links(staad_property=staad_property, geometry=geometry, nodes=nodes, top_nodes_cylinder=top_nodes, nodes_w_constraints=nodes_w_constraints, prop_num=prop_num)
    create_top_contraints(staad_property=staad_property,hole_boundaries=top_nodes, master=nodes_w_constraints[-1])


    load = openstaad.Load
    load._FlagAsMethod("SetLoadActive")
    load._FlagAsMethod("CreateNewPrimaryLoad")
    load._FlagAsMethod("AddSelfWeightInXYZ")
    load._FlagAsMethod("AddNodalLoad")
    case_num = load.CreateNewPrimaryLoad("Loads")
    ret = load.SetLoadActive(case_num)  # Load Case 1
    create_point_loads(load=load,case_num=case_num,load_list=load_list,top_cylinder_nodetag=nodes_w_constraints[-1])

    command = openstaad.Command
    command._FlagAsMethod("PerformAnalysis")
    openstaad._FlagAsMethod("SetSilentMode")
    openstaad._FlagAsMethod("Analyze")
    openstaad._FlagAsMethod("isAnalyzing")
    command.PerformAnalysis(6)
    openstaad.SaveModel(1)
    time.sleep(3)
    openstaad.SetSilentMode(1)

    openstaad.Analyze()
    while openstaad.isAnalyzing():
        print("...Analyzing")
        time.sleep(2)

    output = openstaad.Output
    output._FlagAsMethod("GetPlateCenterVonMisesStresses")
    output._FlagAsMethod("GetNodeDisplacements")
    output._FlagAsMethod("GetSupportReactions")

    von_mises_results = get_all_von_misses(output=output, triangles=triangles, lc=case_num)

    disp_results = get_all_displacements(output=output, triangles=triangles, lc=case_num)
    bolt_reactions = get_all_reactions(output=output, nodes_w_constraints=nodes_w_constraints,lc=1)
    json_path = Path.cwd() / "output.json"
    with open(json_path, "w") as jsonfile:
        json.dump({"von_misses": von_mises_results, "displacement": disp_results, "bolt_reactions":bolt_reactions}, jsonfile)

    openstaad = None
    CoUninitialize()
    staad_process.terminate()
    return None


if __name__ == "__main__":
    run_staad()
