import random 
import viktor as vkt

def generate_colors():
    return random.choice(
        [
            vkt.Color(r=255, g=255, b=0),  # Yellow (Contrasts well with blue/purple)
            vkt.Color(r=0, g=0, b=255),    # Blue (Contrasts well with yellow/orange)
            vkt.Color(r=255, g=69, b=0),   # Red-Orange (Contrasts with cyan/blue)
            vkt.Color(r=0, g=255, b=255),  # Cyan (Contrasts well with red/orange)
            vkt.Color(r=255, g=0, b=255),  # Magenta (Contrasts with green)
        ]
    )


def get_stress_color(stress_value: float) -> vkt.Color:
    """
    Returns a vkt.Color object corresponding to the given stress_value using a custom Jet colormap.
    No interpolation is used—each value falls into a fixed category.
    """
    if stress_value <= 0.190955:
        return vkt.Color(r=0, g=0, b=128)   # Dark blue
    elif stress_value <= 5:
        return vkt.Color(r=0, g=0, b=255)   # Blue
    elif stress_value <= 10:
        return vkt.Color(r=0, g=128, b=255) # Light blue
    elif stress_value <= 15:
        return vkt.Color(r=0, g=255, b=255) # Cyan
    elif stress_value <= 20:
        return vkt.Color(r=0, g=255, b=128) # Greenish cyan
    elif stress_value <= 25:
        return vkt.Color(r=0, g=255, b=0)   # Green
    elif stress_value <= 30:
        return vkt.Color(r=128, g=255, b=0) # Yellow-green
    elif stress_value <= 35:
        return vkt.Color(r=255, g=255, b=0)  # Yellow
    elif stress_value <= 45:
        return vkt.Color(r=255, g=192, b=0)  # Orange-yellow
    elif stress_value <= 55:
        return vkt.Color(r=255, g=128, b=0)  # Orange
    elif stress_value <= 65:
        return vkt.Color(r=255, g=64, b=0)  # Reddish-orange
    elif stress_value <= 90:
        return vkt.Color(r=255, g=0, b=0)  # Red
    elif stress_value <= 120:
        return vkt.Color(r=192, g=0, b=0)  # Dark red
    elif stress_value <= 180:
        return vkt.Color(r=128, g=0, b=0)  # Deep red
    elif stress_value <= 200:
        return vkt.Color(r=128, g=0, b=64)  # Dark magenta
    elif stress_value <= 250:
        return vkt.Color(r=128, g=0, b=128)  # Purple
    else:
        return vkt.Color(r=192, g=0, b=192)  # Pinkish purple

def render_model(nodes: dict, triangles: dict, mesh_mode: bool)->vkt.Group:
        triangle_assemblies = []
        groups = {"cylinder": [], "stiffeners": [], "base plate": []}
        colors = {"cylinder": vkt.Color(r=0, g=255, b=255), "stiffeners": vkt.Color(r=212, g=175, b=55), "base plate": vkt.Color(r=184, g=115, b=51)}

        for triangle_id, triangle_args in triangles.items():
            node_i_tag = triangle_args["NodeI"]
            node_j_tag = triangle_args["NodeJ"]
            node_k_tag = triangle_args["NodeK"]

            node_j = nodes[str(node_j_tag)]
            node_i = nodes[str(node_i_tag)]
            node_k = nodes[str(node_k_tag)]

            point_i = vkt.Point(x=node_i["x"] / 100, y=node_i["y"] / 100, z=node_i["z"] / 100)
            point_j = vkt.Point(x=node_j["x"] / 100, y=node_j["y"] / 100, z=node_j["z"] / 100)
            point_k = vkt.Point(x=node_k["x"] / 100, y=node_k["y"] / 100, z=node_k["z"] / 100)

            triangle_1 = vkt.Triangle(point1=point_i, point2=point_j, point3=point_k)
            triangle_2 = vkt.Triangle(point1=point_k, point2=point_j, point3=point_i)

            if mesh_mode:
                triangle_assemblies.append(
                    vkt.TriangleAssembly(triangles=[triangle_1, triangle_2], skip_duplicate_vertices_check=True, material=vkt.Material(color=generate_colors()))
                )
            else:
                triangle_group = triangle_args.get("group")
                if triangle_group in groups:
                    groups[triangle_group].append(triangle_1)
                    groups[triangle_group].append(triangle_2)

        if mesh_mode:
            return triangle_assemblies
        else:
            objects = [
                vkt.TriangleAssembly(triangles=triangle_list, identifier=group_name, skip_duplicate_vertices_check=True, material=vkt.Material(color=colors[group_name]))
                for group_name, triangle_list in groups.items()
            ]
            return objects

def render_model_results(nodes: dict, triangles: dict, von_misses: dict, displacement: dict, factor: float):
    triangle_assemblies = []
    for triangle_id, triangle_args in triangles.items():
        node_i_tag = triangle_args["NodeI"]
        node_j_tag = triangle_args["NodeJ"]
        node_k_tag = triangle_args["NodeK"]

        node_i = nodes[str(node_i_tag)]
        node_j = nodes[str(node_j_tag)]
        node_k = nodes[str(node_k_tag)]

        disp_i = displacement[str(node_i_tag)]
        disp_j = displacement[str(node_j_tag)]
        disp_k = displacement[str(node_k_tag)]

        # Applying displacement * factor to each coordinate
        point_i = vkt.Point(x=(node_i["x"] + disp_i["x"] * factor) / 100, y=(node_i["y"] + disp_i["z"] * factor) / 100, z=(node_i["z"] + disp_i["y"] * factor) / 100)

        point_j = vkt.Point(x=(node_j["x"] + disp_j["x"] * factor) / 100, y=(node_j["y"] + disp_j["z"] * factor) / 100, z=(node_j["z"] + disp_j["y"] * factor) / 100)

        point_k = vkt.Point(x=(node_k["x"] + disp_k["x"] * factor) / 100, y=(node_k["y"] + disp_k["z"] * factor) / 100, z=(node_k["z"] + disp_k["y"] * factor) / 100)

        triangle_1 = vkt.Triangle(point1=point_i, point2=point_j, point3=point_k)
        triangle_2 = vkt.Triangle(point1=point_k, point2=point_j, point3=point_i)

        triangle_assemblies.append(
            vkt.TriangleAssembly(
                triangles=[triangle_1, triangle_2], skip_duplicate_vertices_check=True, material=vkt.Material(color=get_stress_color(von_misses[str(triangle_id)]))
            )
        )
    return triangle_assemblies


def generate_bolts(nodes, nodes_w_constraints,objects,hole_diameter):
    for bolt_id in nodes_w_constraints[:-1]:
        bolt_node = nodes[str(bolt_id)]

        top_point = vkt.Point(x=bolt_node["x"] / 100, y=bolt_node["y"] / 100, z=bolt_node["z"] / 100)
        bottom_point = vkt.Point(x=bolt_node["x"] / 100, y=bolt_node["y"] / 100, z=bolt_node["z"]-200/ 100)
        bolt_line = vkt.Line(top_point, bottom_point)

        objects.append(vkt.CircularExtrusion(diameter=(hole_diameter - 2) / 100,line=bolt_line, material=vkt.Material(color=vkt.Color(r=255, g=140, b=0))))

        washer_point = vkt.Point(x=bolt_node["x"] / 100, y=bolt_node["y"] / 100, z=bolt_node["z"]/ 100 + 20 / 100)
        washer = vkt.Line(top_point, washer_point)
        objects.append(vkt.CircularExtrusion(diameter=(hole_diameter * 2) / 100, line=washer, material=vkt.Material(color=vkt.Color(r=255, g=140, b=0))))

        top_nut_point = vkt.Point(x=bolt_node["x"] / 100, y=bolt_node["y"] / 100, z= bolt_node["z"]/ 100 + 40 / 100)
        top_nut = vkt.Line(washer_point, top_nut_point)
        objects.append(vkt.CircularExtrusion(diameter=(hole_diameter - 2) / 100, line=top_nut, material=vkt.Material(color=vkt.Color(r=255, g=140, b=0))))
    
    return objects


def generate_foundation_block(widht:float, height:float, depth:float):
    return vkt.RectangularExtrusion(
                    width=widht / 100,
                    height=height / 100,
                    line=vkt.Line(vkt.Point(0, 0, -20 / 100), vkt.Point(0, 0, -depth / 100)),
                    material=vkt.Material(color=vkt.Color(r=128, g=128, b=128), opacity=0.7),
                )
