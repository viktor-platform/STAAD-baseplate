import math
import viktor as vkt
from munch import Munch


load_dir2vec = {
    "Y+": vkt.Vector(0, 1, 0),
    "Y-": vkt.Vector(0, -1, 0),
    "X-": vkt.Vector(-1, 0, 0),
    "X+": vkt.Vector(1, 0, 0),
    "Z+": vkt.Vector(0, 0, 1),
    "Z-": vkt.Vector(0, 0, -1)
    }


def process_loads(load_list:list[Munch]) -> list[dict]:
    """
    Function that process a list of loads by creating a dict that will be sent to staad.
    """
    loads_dict = []
    for load in load_list:
        mag = load.magnitud
        ld_dir = load.dir
        load_type = load.type

        loads_dict.append({"type": load_type, "mag": mag, "ld_dir": ld_dir,"color":load.color})
    return loads_dict


def create_arrow_at_point(tip: vkt.Point, direction: vkt.Vector, arrow_size: float, material: vkt.Material)->vkt.Group:
    """
    Function taht creates an arrow (shaft and head) for vertical and horizontal loads so that
    the pointed tip of the arrow head touches the given point.
    """
    arrow_head_height = arrow_size / 2.0
    arrow_head_radius = arrow_size / 4.0

    # The shaft extends from a tail (offset by the full arrow_size)
    # to the base of the arrow head (offset by arrow_head_height)
    tail = vkt.Point(
        tip.x - direction.x * arrow_size,
        tip.y - direction.y * arrow_size,
        tip.z - direction.z * arrow_size,
    )
    arrow_head_base = vkt.Point(
        tip.x - direction.x * arrow_head_height,
        tip.y - direction.y * arrow_head_height,
        tip.z - direction.z * arrow_head_height,
    )

    shaft_thickness = arrow_size / 10.0
    arrow_shaft = vkt.RectangularExtrusion(shaft_thickness, shaft_thickness, vkt.Line(tail, arrow_head_base), material=material)

    arrow_head = vkt.Cone(arrow_head_radius, arrow_head_height, origin=arrow_head_base, orientation=direction, material=material)

    return vkt.Group([arrow_shaft, arrow_head])


def create_moment_arrow_head(arc_end: vkt.Point, direction: vkt.Vector, arrow_size: float, material: vkt.Material)->vkt.Cone:
    """Function to create the moment arrow head"""
    arrow_head_height = arrow_size / 2.0
    arrow_head_radius = arrow_size / 4.0
    # Use the arc end point as the cone's origin.
    return vkt.Cone(arrow_head_radius, arrow_head_height, origin=arc_end, orientation=direction, material=material)


def create_moment_loads(x: float, y: float, z: float, arrow_size: float, color: vkt.Color)->vkt.Group:
    """ Function that create the moment arrow to be render in the GeometryView"""
    moment_elements = []
    radius = 200 / 100
    # Define a red material.
    material = vkt.Material(color=color)
    num_segments = 20
    arc_points = []
    # Generate the upper half circle in the x–z plane.
    for i in range(num_segments + 1):
        theta = i * math.pi / num_segments  # theta from 0 to pi
        pt = vkt.Point(x + radius * math.cos(theta), y, z + radius * math.sin(theta))
        arc_points.append(pt)

    extruded_segments = []
    arc_thickness = arrow_size / 10.0
    for i in range(len(arc_points) - 1):
        segment_line = vkt.Line(arc_points[i], arc_points[i + 1])
        extruded_segment = vkt.RectangularExtrusion(arc_thickness, arc_thickness, segment_line, material=material)
        extruded_segments.append(extruded_segment)

    moment_arc = vkt.Group(extruded_segments)
    moment_elements.append(moment_arc)
    arc_end = arc_points[0]
    moment_arrow_direction = vkt.Vector(0, 0, -1)  # Arrow head points downward.
    moment_arrow_head = create_moment_arrow_head(arc_end, moment_arrow_direction, arrow_size, material)
    moment_elements.append(moment_arrow_head)
    return moment_elements


def create_point_loads(x: float, y: float, z:float , arrow_size: float, dir_vector: vkt.Vector ,color: vkt.Color)->vkt.Group:
    """Function that creates the point loads to be render in the geometry view"""
    material = vkt.Material(color=color)
    offset = 10/100
    node = vkt.Point(x, y, z + offset)
    vertical_arrow = create_arrow_at_point(node, dir_vector, arrow_size, material)
    return vertical_arrow
