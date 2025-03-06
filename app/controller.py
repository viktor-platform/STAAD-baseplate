import json
import math
import viktor as vkt

from app.baseplate import TubularBasePlateModel
from app.parametrization import Parametrization
from app.visualization import generate_bolts, generate_foundation_block, render_model, render_model_results
from app.loads import create_moment_loads, create_point_loads, load_dir2vec, process_loads

from io import BytesIO
from pathlib import Path
from viktor.core import File
from viktor.external.python import PythonAnalysis


def get_model_args(params):
        return  {
            "cylinder_radius": params.step_1.cylinder_radius,
            "cylinder_height": params.step_1.cylinder_height,
            "num_circ_div": params.step_1.num_circ_div,
            "num_height_div": params.step_1.num_height_div,
            "num_stiffeners": params.step_1.num_stiffeners,
            "stiffener_height": params.step_1.stiffener_height,
            "stiffener_base": params.step_1.stiffener_base,
            "num_stiffener_divisions": params.step_1.num_stiffener_divisions,
            "baseplate_radius": params.step_1.baseplate_radius,
            "num_baseplate_div": params.step_1.num_baseplate_div,
            "num_holes": params.step_1.num_holes,
            "hole_margin": params.step_1.hole_margin,
            "hole_angle_offset": math.pi / params.step_1.num_stiffeners,
            "hole_diameter": params.step_1.hole_diameter,
        }

@vkt.memoize
def generate_model(kwargs) -> tuple[dict[str, dict], dict[str, dict], list[str]]:
    model = TubularBasePlateModel(**kwargs)
    nodes, triangles, supports = model.build()
    hole_boundaries = model.get_hole_boundary_node_tags()
    top_nodes = model.get_top_node_tags()
    return nodes, triangles, supports, hole_boundaries, top_nodes


@vkt.memoize
def run_staad_analysis(kwargs, plate_thickness:float, loads: dict|None = None):
    nodes, triangles, supports, hole_boundaries, top_nodes = generate_model(kwargs)
    # Remove vkt.color for serialization
    if loads:
        [d.pop("color", None) for d in loads]
    else:
        loads = []

    input_json = json.dumps([nodes, triangles, supports, hole_boundaries, top_nodes,loads,plate_thickness])
    script_path = Path(__file__).parent / "run_staad_model.py"
    macro_path = Path(__file__).parent / "thickness.xlsm"
    script = File.from_path(script_path)

    files = [("inputs.json", BytesIO(bytes(input_json, "utf8"))), ("thickness.xlsm", File.from_path(macro_path))]

    staad_analysis = PythonAnalysis(script=script, files=files, output_filenames=["output.json"])
    staad_analysis.execute(timeout=3600)
    output_file = staad_analysis.get_output_file("output.json")
    output_file = json.loads(output_file.getvalue())
    return output_file["von_misses"], output_file["displacement"], output_file["bolt_reactions"]


class Controller(vkt.Controller):
    parametrization = Parametrization(width=30)

    @vkt.GeometryView(label="Base Plate 3D Model", duration_guess=40, x_axis_to_right=True)
    def render_model(self, params, **kwargs):
        
        model_kwargs = get_model_args(params)
        nodes, triangles, nodes_w_constraints, hole_boundaries, top_nodes = generate_model(model_kwargs)

        objects = render_model(nodes=nodes,triangles=triangles,mesh_mode=params.step_1.mode)

        # Render bolts and foundation block
        if not params.step_1.mode:
            objects = generate_bolts(nodes=nodes, nodes_w_constraints=nodes_w_constraints, objects=objects, hole_diameter=params.step_1.hole_diameter)

            objects.append(
                generate_foundation_block(
                    widht=params.step_1.baseplate_radius*2*1.1,
                    height=params.step_1.baseplate_radius*2*1.1,
                    depth= params.step_1.baseplate_radius*2*1.5,
                ))

        # Render loads
        load_point = nodes_w_constraints[-1]
        lp_coords = nodes[str(load_point)]
        loads = params.step_2.loads if params.step_2.loads else None
        if loads:
            load_list = process_loads(loads)
            for load_dict in load_list:
                if load_dict["type"] == "Moment Load":
                    moment_elements = create_moment_loads(lp_coords["x"] / 100, lp_coords["y"] / 100, lp_coords["z"] / 100, arrow_size=1,color=load_dict["color"])
                    objects.extend(moment_elements)
                else:
                    point_load = create_point_loads(lp_coords["x"] / 100, lp_coords["y"] / 100, lp_coords["z"] / 100, arrow_size=1, dir_vector=load_dir2vec[load_dict["ld_dir"]],color=load_dict["color"])
                    objects.append(point_load)

        return vkt.GeometryResult(geometry=objects)
        
    @vkt.GeometryView(label="Stress Distributions", duration_guess=40, x_axis_to_right=True)
    def results(self, params, **kwargs):
        model_kwargs = get_model_args(params)
        loads = params.step_2.loads if params.step_2.loads else None
        if loads:
            load_list = process_loads(loads)
        else:
            load_list = []
        von_misses, displacement, reactions = run_staad_analysis(kwargs=model_kwargs,plate_thickness=params.step_1.plate_thickness,loads=load_list)
        nodes, triangles, supports, hole_boundaries, top_nodes = generate_model(model_kwargs)
        factor = params.step_3.scale
        triangle_assemblies = render_model_results(nodes=nodes, triangles=triangles, von_misses=von_misses, displacement=displacement, factor=factor)
        return vkt.GeometryResult(geometry=triangle_assemblies)


