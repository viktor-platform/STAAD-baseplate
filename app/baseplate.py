import math
import numpy as np
import meshpy.triangle as triangle

class Cylinder:
    @staticmethod
    def create_mesh(
        radius: float,
        height: float,
        num_theta_pts: int,
        num_z_pts: int,
        extra_theta: list[float] | None = None,
        extra_z: list[float] | None = None
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        tol: float = 1e-8
        regular_theta: list[float] = [2 * math.pi * i / num_theta_pts for i in range(num_theta_pts)]
        theta_set: list[float] = list(regular_theta)
        if extra_theta is not None:
            for a in extra_theta:
                a_mod = a % (2 * math.pi)
                if not any(abs(a_mod - t) < tol for t in theta_set):
                    theta_set.append(a_mod)
        theta_set.sort()

        regular_z: list[float] = np.linspace(0, height, num_z_pts + 1).tolist()
        z_set: list[float] = list(regular_z)
        if extra_z is not None:
            for z in extra_z:
                if 0 <= z <= height and not any(abs(z - z0) < tol for z0 in z_set):
                    z_set.append(z)
        z_set.sort()

        n_theta: int = len(theta_set)
        n_z: int = len(z_set)
        pts: list[list[float]] = []
        params: list[list[float]] = []
        for th in theta_set:
            for z in z_set:
                x: float = radius * math.cos(th)
                y: float = radius * math.sin(th)
                pts.append([x, y, z])
                params.append([th, z])
        pts_array: np.ndarray = np.array(pts)

        elems: list[tuple[int, int, int]] = []
        for i in range(n_theta):
            next_i: int = (i + 1) % n_theta  # periodic
            for j in range(n_z - 1):
                idx0: int = i * n_z + j
                idx1: int = next_i * n_z + j
                idx2: int = i * n_z + (j + 1)
                idx3: int = next_i * n_z + (j + 1)
                elems.append((idx0, idx1, idx2))
                elems.append((idx1, idx3, idx2))
        elems_array: np.ndarray = np.array(elems)
        return pts_array, elems_array, np.array(params)

    def __init__(
        self,
        radius: float,
        height: float,
        num_circ_div: int,
        num_height_div: int,
        extra_theta: list[float] | None = None,
        extra_z: list[float] | None = None
    ) -> None:
        self.radius = radius
        self.height = height
        self.num_circ_div = num_circ_div
        self.num_height_div = num_height_div
        self.extra_theta = extra_theta
        self.extra_z = extra_z

    def generate_mesh(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return Cylinder.create_mesh(
            self.radius,
            self.height,
            self.num_circ_div,
            self.num_height_div,
            self.extra_theta,
            self.extra_z
        )

    def get_top_node_ids(self, mesh_nodes: np.ndarray, tol: float = 1e-6) -> list[int]:
        """Return indices of nodes with z-coordinate equal (within tol) to the cylinder height."""
        top_ids: list[int] = []
        for i, node in enumerate(mesh_nodes):
            if abs(node[2] - self.height) < tol:
                top_ids.append(i)
        return top_ids


class Stiffener:
    @staticmethod
    def create_mesh(
        angle: float,
        cylinder_radius: float,
        height: float,
        base_length: float,
        num_divisions: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        N: int = num_divisions
        pts2d: list[tuple[float, float]] = []
        index_map: dict[tuple[int, int], int] = {}
        count: int = 0
        for j in range(N + 1):
            for i in range(N - j + 1):
                px: float = base_length * (i / N)
                py: float = height * (j / N)
                pts2d.append((px, py))
                index_map[(i, j)] = count
                count += 1
        pts2d_array: np.ndarray = np.array(pts2d)

        elems: list[tuple[int, int, int]] = []
        for j in range(N):
            n_in_row: int = (N - j + 1)
            for i in range(n_in_row - 1):
                a: int = index_map[(i, j)]
                b: int = index_map[(i + 1, j)]
                c: int = index_map[(i, j + 1)]
                elems.append((a, b, c))
                if i < (N - j - 1):
                    d: int = index_map[(i + 1, j + 1)]
                    elems.append((b, d, c))
        elems_array: np.ndarray = np.array(elems)

        pts3d: list[list[float]] = []
        for (px, py) in pts2d:
            base_x: float = cylinder_radius * math.cos(angle)
            base_y: float = cylinder_radius * math.sin(angle)
            offset_x: float = px * math.cos(angle)
            offset_y: float = px * math.sin(angle)
            pts3d.append([base_x + offset_x, base_y + offset_y, py])
        pts3d_array: np.ndarray = np.array(pts3d)
        return pts3d_array, elems_array, pts2d_array

    @staticmethod
    def snap_to_cylinder(
        stiff_3d: np.ndarray,
        stiff_elem: np.ndarray,
        stiff_param: np.ndarray,
        angle: float,
        cylinder_radius: float,
        param_tol: float = 1e-6
    ) -> tuple[np.ndarray, np.ndarray]:
        stiff_3d_out: np.ndarray = stiff_3d.copy()
        for i, (px, py) in enumerate(stiff_param):
            if abs(px) < param_tol:
                stiff_3d_out[i] = [cylinder_radius * math.cos(angle),
                                   cylinder_radius * math.sin(angle),
                                   py]
        return stiff_3d_out, stiff_elem

    @staticmethod
    def snap_to_baseplate(
        stiff_3d: np.ndarray,
        stiff_param: np.ndarray,
        base_length: float,
        baseplate_radius: float,
        angle: float,
        tol: float = 1e-6
    ) -> np.ndarray:
        stiff_3d_out: np.ndarray = stiff_3d.copy()
        for i, (px, py) in enumerate(stiff_param):
            if abs(py) < tol and abs(px - base_length) < tol:
                stiff_3d_out[i] = [baseplate_radius * math.cos(angle),
                                   baseplate_radius * math.sin(angle),
                                   0.0]
        return stiff_3d_out

    def __init__(
        self,
        angle: float,
        cylinder_radius: float,
        stiffener_height: float,
        stiffener_base: float,
        num_stiffener_divisions: int
    ) -> None:
        self.angle = angle
        self.cylinder_radius = cylinder_radius
        self.stiffener_height = stiffener_height
        self.stiffener_base = stiffener_base
        self.num_stiffener_divisions = num_stiffener_divisions

    def generate_mesh(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        stf_3d, stf_elem, stf_param = Stiffener.create_mesh(
            self.angle,
            self.cylinder_radius,
            self.stiffener_height,
            self.stiffener_base,
            self.num_stiffener_divisions
        )
        stf_3d, stf_elem = Stiffener.snap_to_cylinder(
            stf_3d, stf_elem, stf_param, self.angle, self.cylinder_radius
        )
        stf_3d = Stiffener.snap_to_baseplate(
            stf_3d, stf_param, self.stiffener_base, self.cylinder_radius + self.stiffener_base, self.angle
        )
        return stf_3d, stf_elem, stf_param


class CylinderWithStiffeners(Cylinder):
    @staticmethod
    def merge_meshes(
        base_nodes: np.ndarray,
        base_elems: list[tuple[int, int, int]],
        base_groups: list[str],
        add_nodes: np.ndarray,
        add_elems: np.ndarray,
        add_groups: list[str],
        tol: float = 1e-8
    ) -> tuple[np.ndarray, np.ndarray, list[str]]:
        merged_nodes: list[np.ndarray] = list(base_nodes)
        index_map: dict[int, int] = {}
        for i, n in enumerate(add_nodes):
            found_match: bool = False
            for j, m in enumerate(merged_nodes):
                if np.linalg.norm(n - m) < tol:
                    index_map[i] = j
                    found_match = True
                    break
            if not found_match:
                index_map[i] = len(merged_nodes)
                merged_nodes.append(n)
        merged_nodes_array: np.ndarray = np.array(merged_nodes)
        merged_add_elems: list[tuple[int, int, int]] = []
        for tri in add_elems:
            ia: int = index_map[tri[0]]
            ib: int = index_map[tri[1]]
            ic: int = index_map[tri[2]]
            merged_add_elems.append((ia, ib, ic))
        merged_elems: list[tuple[int, int, int]] = list(base_elems) + merged_add_elems
        merged_groups: list[str] = base_groups + add_groups
        return merged_nodes_array, np.array(merged_elems), merged_groups

    def __init__(self, cylinder: Cylinder, stiffeners: list[Stiffener]) -> None:
        self.__dict__.update(cylinder.__dict__)
        self.stiffeners = stiffeners

    def generate_mesh(self) -> tuple[np.ndarray, np.ndarray, list[str]]:
        cyl_nodes, cyl_elems, _ = super().generate_mesh()
        cyl_groups: list[str] = ["cylinder"] * len(cyl_elems)
        merged_nodes: np.ndarray = cyl_nodes
        merged_elems: list[tuple[int, int, int]] = list(cyl_elems)
        merged_groups: list[str] = list(cyl_groups)
        for stiff in self.stiffeners:
            stf_nodes, stf_elems, _ = stiff.generate_mesh()
            groups: list[str] = ["stiffeners"] * len(stf_elems)  # Ensure group label is "stiffeners"
            merged_nodes, merged_elems, merged_groups = CylinderWithStiffeners.merge_meshes(
                merged_nodes, merged_elems, merged_groups,
                stf_nodes, stf_elems, groups,
                tol=1e-8
            )
        return merged_nodes, np.array(merged_elems), merged_groups


class CircularPlate:
    @staticmethod
    def generate_circle_polygon(
        center: tuple[float, float],
        diameter: float,
        num_points: int = 20
    ) -> list[tuple[float, float]]:
        cx, cy = center
        r: float = diameter / 2.0
        pts: list[tuple[float, float]] = []
        for i in range(num_points):
            theta: float = 2 * math.pi * i / num_points
            pts.append((cx + r * math.cos(theta), cy + r * math.sin(theta)))
        return pts

    @staticmethod
    def generate_radial_hole_centers(
        num_holes: int,
        baseplate_radius: float,
        hole_margin: float,
        angle_offset: float = 0.0
    ) -> list[tuple[float, float]]:
        centers: list[tuple[float, float]] = []
        if num_holes == 1:
            centers.append((0.0, 0.0))
        else:
            r: float = baseplate_radius - hole_margin
            for i in range(num_holes):
                angle: float = angle_offset + 2 * math.pi * i / num_holes
                centers.append((r * math.cos(angle), r * math.sin(angle)))
        return centers

    @staticmethod
    def create_mesh(
        baseplate_radius: float,
        num_div: int,
        num_holes: int,
        hole_diameter: float = 12,
        mesh_max_area: float | None = None,
        extra_points: list[tuple[float, float]] | None = None,
        hole_margin: float = 10,
        hole_angle_offset: float = 0.0
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        outer: list[tuple[float, float]] = CircularPlate.generate_circle_polygon((0.0, 0.0), diameter=2 * baseplate_radius, num_points=40)
        pts: list[tuple[float, float]] = outer.copy()
        segments: list[tuple[int, int]] = [(i, (i + 1) % len(outer)) for i in range(len(outer))]
        hole_centers: list[tuple[float, float]] = CircularPlate.generate_radial_hole_centers(num_holes, baseplate_radius, hole_margin, angle_offset=hole_angle_offset)
        for hc in hole_centers:
            poly: list[tuple[float, float]] = CircularPlate.generate_circle_polygon(hc, hole_diameter, num_points=20)
            start_idx: int = len(pts)
            pts.extend(poly)
            num_poly: int = len(poly)
            segs: list[tuple[int, int]] = [(start_idx + i, start_idx + ((i + 1) % num_poly)) for i in range(num_poly)]
            segments.extend(segs)
        if extra_points is not None:
            pts.extend(extra_points)
        mesh_info: triangle.MeshInfo = triangle.MeshInfo()
        mesh_info.set_points(pts)
        mesh_info.set_facets(segments)
        mesh_info.set_holes(hole_centers)
        target_area: float = mesh_max_area if mesh_max_area is not None else (math.pi * (baseplate_radius ** 2)) / (num_div ** 2)
        def refine_func(vertices: list[tuple[float, float]], area: float) -> bool:
            return area > target_area
        mesh2d: triangle.MeshInfo = triangle.build(mesh_info, refinement_func=refine_func)
        pts3d: list[list[float]] = []
        for x, y in mesh2d.points:
            pts3d.append([x, y, 0.0])
        return np.array(pts3d), np.array(mesh2d.elements), np.array(mesh2d.points)

    def __init__(
        self,
        baseplate_radius: float,
        num_baseplate_div: int,
        num_holes: int,
        hole_diameter: float,
        hole_margin: float,
        hole_angle_offset: float
    ) -> None:
        self.baseplate_radius = baseplate_radius
        self.num_baseplate_div = num_baseplate_div
        self.num_holes = num_holes
        self.hole_diameter = hole_diameter
        self.hole_margin = hole_margin
        self.hole_angle_offset = hole_angle_offset

    def generate_mesh(self, extra_points: list[tuple[float, float]] | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return CircularPlate.create_mesh(
            self.baseplate_radius,
            self.num_baseplate_div,
            self.num_holes,
            hole_diameter=self.hole_diameter,
            mesh_max_area=None,
            extra_points=extra_points,
            hole_margin=self.hole_margin,
            hole_angle_offset=self.hole_angle_offset
        )

    def get_hole_boundary_node_ids(self, mesh_nodes: np.ndarray, tol: float = 1e-3) -> dict[str, list[int]]:
        """
        Returns a dictionary where each key is "hole{i}" and its value is a list of node indices (from mesh_nodes)
        that lie on the perimeter of that hole.
        """
        boundaries: dict[str, list[int]] = {}
        expected_radius: float = self.hole_diameter / 2.0
        hole_centers: list[tuple[float, float]] = CircularPlate.generate_radial_hole_centers(
            self.num_holes, self.baseplate_radius, self.hole_margin, self.hole_angle_offset
        )
        for i, center in enumerate(hole_centers):
            ids: list[int] = []
            for j, node in enumerate(mesh_nodes):
                if abs(node[2]) < 1e-6:  # Only consider baseplate nodes
                    dist = math.hypot(node[0] - center[0], node[1] - center[1])
                    if abs(dist - expected_radius) < tol:
                        ids.append(j)
            boundaries[f"hole{i}"] = ids
        return boundaries


class BaseModel:
    @staticmethod
    def serialize_mesh(
        nodes: np.ndarray,
        triangles: np.ndarray,
        tri_groups: list[str],
        filename: str = "mesh.json",
        extra_node_ids: list[int] | None = None
    ) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]], list[int] | None]:
        mesh_data: dict[str, dict[str, dict[str, float]]] = {"nodes": {}, "triangles": {}}
        for i, node in enumerate(nodes):
            mesh_data["nodes"][str(i + 1)] = {
                "id": i + 1,
                "x": node[0],
                "y": node[1],
                "z": node[2]
            }
        for i, tri in enumerate(triangles):
            mesh_data["triangles"][str(i + 1)] = {
                "id": i + 1,
                "NodeI": int(tri[0]) + 1,
                "NodeJ": int(tri[1]) + 1,
                "NodeK": int(tri[2]) + 1,
                "group": tri_groups[i]
            }
        if extra_node_ids is not None:
            mesh_data["extra_node_ids"] = [nid + 1 for nid in extra_node_ids]
        return mesh_data["nodes"], mesh_data["triangles"], mesh_data.get("extra_node_ids", None)

    @staticmethod
    def clean_model(
        Nodes: dict[str, dict[str, float]],
        Triangles: dict[str, dict[str, float]]
    ) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]]]:
        tol: float = 1e-8
        seen: dict[tuple[float, float, float], int] = {}
        new_nodes: dict[int, dict[str, float]] = {}
        mapping: dict[int, int] = {}
        new_id: int = 1
        for key in sorted(Nodes.keys(), key=lambda k: int(k)):
            node: dict[str, float] = Nodes[key]
            coord: tuple[float, float, float] = (round(node["x"], 8), round(node["y"], 8), round(node["z"], 8))
            if coord in seen:
                mapping[node["id"]] = seen[coord]
            else:
                seen[coord] = new_id
                mapping[node["id"]] = new_id
                new_nodes[new_id] = {
                    "id": new_id,
                    "x": node["x"],
                    "y": node["y"],
                    "z": node["z"]
                }
                new_id += 1
        new_triangles: dict[int, dict[str, float]] = {}
        new_tri_id: int = 1
        for key in sorted(Triangles.keys(), key=lambda k: int(k)):
            tri: dict[str, float] = Triangles[key]
            ni: int = mapping[tri["NodeI"]]
            nj: int = mapping[tri["NodeJ"]]
            nk: int = mapping[tri["NodeK"]]
            new_triangles[new_tri_id] = {
                "id": new_tri_id,
                "NodeI": ni,
                "NodeJ": nj,
                "NodeK": nk,
                "group": tri["group"]
            }
            new_tri_id += 1
        cleaned_nodes: dict[str, dict[str, float]] = {str(k): v for k, v in new_nodes.items()}
        cleaned_triangles: dict[str, dict[str, float]] = {str(k): v for k, v in new_triangles.items()}
        return cleaned_nodes, cleaned_triangles


class TubularBasePlateModel(BaseModel):
    def __init__(
        self,
        cylinder_radius: float,
        cylinder_height: float,
        num_circ_div: int,
        num_height_div: int,
        num_stiffeners: int,
        stiffener_height: float,
        stiffener_base: float,
        num_stiffener_divisions: int,
        baseplate_radius: float,
        num_baseplate_div: int,
        num_holes: int,
        hole_diameter: float,
        hole_margin: float,
        hole_angle_offset: float
    ) -> None:
        extra_theta: list[float] = [2 * math.pi * i / num_stiffeners for i in range(num_stiffeners)]
        extra_z: list[float] = [j * stiffener_height / num_stiffener_divisions for j in range(num_stiffener_divisions + 1)]
        self.cylinder = Cylinder(
            radius=cylinder_radius,
            height=cylinder_height,
            num_circ_div=num_circ_div,
            num_height_div=num_height_div,
            extra_theta=extra_theta,
            extra_z=extra_z
        )
        self.stiffeners: list[Stiffener] = []
        for i in range(num_stiffeners):
            angle: float = 2 * math.pi * i / num_stiffeners if num_stiffeners > 1 else 0.0
            self.stiffeners.append(
                Stiffener(
                    angle=angle,
                    cylinder_radius=cylinder_radius,
                    stiffener_height=stiffener_height,
                    stiffener_base=stiffener_base,
                    num_stiffener_divisions=num_stiffener_divisions
                )
            )
        self.cyl_with_stiffs = CylinderWithStiffeners(self.cylinder, self.stiffeners)
        self.circular_plate = CircularPlate(
            baseplate_radius=baseplate_radius,
            num_baseplate_div=num_baseplate_div,
            num_holes=num_holes,
            hole_diameter=hole_diameter,
            hole_margin=hole_margin,
            hole_angle_offset=hole_angle_offset
        )
        self.cylinder_height = cylinder_height
        self._merged_nodes: np.ndarray | None = None  # Will store merged nodes after build

    def build(self) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]], list[int] | None]:
        # Generate mesh from cylinder with stiffeners.
        cyl_nodes, cyl_elems, cyl_groups = self.cyl_with_stiffs.generate_mesh()
        # Determine forced nodes (those at z ≈ 0) to enforce interface.
        forced_nodes: set[tuple[float, float]] = set()
        for p in cyl_nodes:
            if abs(p[2]) < 1e-6:
                forced_nodes.add((p[0], p[1]))
        # Generate the circular plate mesh using forced nodes.
        baseplate_nodes, baseplate_elems, _ = self.circular_plate.generate_mesh(list(forced_nodes))
        baseplate_groups: list[str] = ["base plate"] * len(baseplate_elems)
        merged_nodes, merged_elems, merged_groups = CylinderWithStiffeners.merge_meshes(
            cyl_nodes, list(cyl_elems), cyl_groups,
            baseplate_nodes, baseplate_elems, baseplate_groups,
            tol=1e-8
        )
        # Append extra support nodes in a fixed order:
        # First: bolt centers (middle of each hole, in order).
        # Last: cylinder top center.
        bolt_centers: list[tuple[float, float]] = CircularPlate.generate_radial_hole_centers(
            self.circular_plate.num_holes,
            self.circular_plate.baseplate_radius,
            self.circular_plate.hole_margin,
            angle_offset=self.circular_plate.hole_angle_offset
        )
        extra_nodes_list: list[np.ndarray] = [np.array([x, y, 0.0]) for (x, y) in bolt_centers]
        cylinder_top_node: np.ndarray = np.array([0.0, 0.0, self.cylinder_height])
        extra_nodes_list.append(cylinder_top_node)
        start_extra: int = merged_nodes.shape[0]
        merged_nodes = np.vstack([merged_nodes, np.array(extra_nodes_list)])
        extra_node_ids: list[int] = list(range(start_extra, merged_nodes.shape[0]))
        nodes, triangles, supports = self.serialize_mesh(
            merged_nodes, merged_elems, merged_groups, filename="mesh.json", extra_node_ids=extra_node_ids
        )
        cleaned_nodes, cleaned_triangles = self.clean_model(nodes, triangles)
        self._merged_nodes = merged_nodes
        return cleaned_nodes, cleaned_triangles, supports

    def get_hole_boundary_node_tags(self, tol: float = 1e-3) -> dict[str, list[int]]:
        """
        Returns a dictionary with keys 'hole0', 'hole1', ... corresponding to each hole,
        and each value is a list of node indices (from the merged nodes) forming that hole's perimeter.
        """
        if self._merged_nodes is None:
            raise ValueError("Model has not been built yet. Call build() first.")
        boundaries: dict[str, list[int]] = {}
        expected_radius: float = self.circular_plate.hole_diameter / 2.0
        hole_centers: list[tuple[float, float]] = CircularPlate.generate_radial_hole_centers(
            self.circular_plate.num_holes,
            self.circular_plate.baseplate_radius,
            self.circular_plate.hole_margin,
            self.circular_plate.hole_angle_offset
        )
        for i, center in enumerate(hole_centers):
            ids: list[int] = []
            for j, node in enumerate(self._merged_nodes):
                if abs(node[2]) < 1e-6:
                    dist = math.hypot(node[0] - center[0], node[1] - center[1])
                    if abs(dist - expected_radius) < tol:
                        ids.append(j)
            boundaries[f"hole{i}"] = ids
        return boundaries

    def get_top_node_tags(self, tol: float = 1e-6) -> dict[str, list[int]]:
        """
        Returns a dictionary with key 'top' whose value is a list of node indices
        from the merged nodes that have z-coordinate equal (within tol) to the cylinder height.
        """
        if self._merged_nodes is None:
            raise ValueError("Model has not been built yet. Call build() first.")
        top_ids: list[int] = []
        for i, node in enumerate(self._merged_nodes):
            if abs(node[2] - self.cylinder_height) < tol:
                top_ids.append(i)
        return {"top": top_ids}