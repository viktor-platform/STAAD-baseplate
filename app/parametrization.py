import viktor as vkt
from textwrap import dedent

suffix = "mm"


class Parametrization(vkt.Parametrization):
    step_1 = vkt.Step(title="Parametric Base Plate - STAAD.PRO", views=["render_model"])
    step_1.text_1 = vkt.Text(
        dedent("""\
            # Create Base Plate Geometry!
            You can change the following input to change the geometry
            and parameters of the base plate as well as the mesh setting!
        """)
    )
    # Cylinder and stiffener parameters.
    step_1.text_2 = vkt.Text(
        dedent("""# Create CHS/Tube Section
    """)
    )
    step_1.cylinder_radius = vkt.NumberField("Cylinder Radius", default=100.0, suffix=suffix)
    step_1.cylinder_height = vkt.NumberField("Cylinder Height", default=300.0, suffix=suffix)
    step_1.lb_1 = vkt.LineBreak()

    step_1.num_circ_div = vkt.NumberField("Radial Mesh Divisions", default=20)
    step_1.num_height_div = vkt.NumberField("Vertical Mesh Divisions", default=20)

    step_1.text_3 = vkt.Text(
        dedent("""# Create Stiffeners
    """)
    )
    step_1.num_stiffeners = vkt.NumberField("Number of Stiffeners", default=7)
    step_1.stiffener_height = vkt.NumberField("Stiffeners Height", default=150.0, suffix=suffix)
    step_1.lb_2 = vkt.LineBreak()

    step_1.stiffener_base = vkt.NumberField("Stiffener Base Length", default=110, suffix=suffix)
    step_1.num_stiffener_divisions = vkt.NumberField("Mesh Divisions", default=10)

    # Circular Baseplate parameters.
    step_1.text_4 = vkt.Text(
        dedent("""# Create Base Plate
    """)
    )
    step_1.baseplate_radius = vkt.NumberField("Base Plate Radius", default=220.0, suffix=suffix)
    step_1.num_baseplate_div = vkt.NumberField("Mesh Divisions", default=30, suffix=suffix)

    step_1.text_5 = vkt.Text(
        dedent("""# Create Base Plate Holes
    """)
    )
    step_1.num_holes = vkt.NumberField("Number of Bolts", default=7)
    step_1.hole_diameter = vkt.NumberField("Bolt Diameter", default=22)
    step_1.lb_3 = vkt.LineBreak()

    step_1.hole_margin = vkt.NumberField("Hole Edge Distance", default=50)
    step_1.text_6 = vkt.Text("# Plate Thickness")
    step_1.plate_thickness = vkt.OptionField("Plate Thickness [mm]", options=[6,10,12,14,16,20],default=12)
    step_1.text_7 = vkt.Text(
        dedent("""# Render Mesh!
    Toggle the component below to display the mesh!
    """)
    )
    step_1.mode = vkt.BooleanField("Display Mesh", default=False)
    step_1.lb_7 = vkt.LineBreak()

    step_2 = vkt.Step(title="Loads", views=["render_model"])
    step_2.text = vkt.Text(dedent("""# Assign the Loads!
    You can create a dynamic array and assign either a moment or
    a point load, defining its direction and magnitude. The loads
    will be rendered in the model.
    """))
    step_2.loads = vkt.DynamicArray("## Add Loads", copylast=True)
    step_2.loads.type = vkt.OptionField("Type", options=["Moment Load", "Point Load"])
    step_2.loads.dir = vkt.OptionField("Direction", options=["Y+", "Y-", "X-", "X+", "Z+", "Z-"])
    step_2.loads.magnitud = vkt.NumberField("Magnitude [kN]", min=0, default=10)
    step_2.loads.color = vkt.ColorField("Color", default=vkt.Color(128, 128, 128))

    step_3 = vkt.Step(title="Results", views=["results"])
    step_3.main_text = vkt.Text("""# Stress Distribution

    The stress distribution is color coded as follows:

    - 🟦 <= 0 [MPa]
    - 🟦 <= 15 [MPa]
    - 🟩 <= 25 [MPa]
    - 🟨 <= 35 [MPa]
    - 🟧 <= 55 [MPa]
    - 🟧 <= 65 [MPa]
    - 🟥 <= 90 [MPa]
    - 🟥 <= 120 [MPa]
    - 🟪 *> 200 [MPa]
    """)
    step_3.br = vkt.LineBreak()
    step_3.text_step2 = vkt.Text("""## Scale Factor
    Move the slider to modify the deformed shape of the 3D model!
    """)
    step_3.scale = vkt.NumberField(""" """, default=100, min=1, max=300, variant="slider")
