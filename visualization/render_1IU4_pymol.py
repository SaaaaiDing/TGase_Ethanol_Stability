"""Ray-traced, aligned 1IU4–ZQG representative figures in PyMOL.

Run in this directory with: pymol -cq render_1IU4_pymol.py
Input PDBs are distributed alongside this script.
"""

from pathlib import Path
from math import dist

from pymol import cmd


HERE = Path(__file__).resolve().parent
OUT = HERE / "pymol_figures"
OUT.mkdir(exist_ok=True)

SYSTEMS = [
    {
        "key": "water", "panel": "A", "title": "Water", "time": "116.0 ns",
        "pdb": "water_116ns.pdb", "residues": (64, 75),
        "contacts": ((64, "N", "O22", 2.73), (75, "OH", "O24", 2.74)),
        "reactive": 3.47,
    },
    {
        "key": "ethanol20", "panel": "B", "title": "20% ethanol", "time": "125.4 ns",
        "pdb": "ethanol20_125p4ns.pdb", "residues": (64, 208, 238),
        "contacts": ((208, "NH2", "O24", 2.77), (238, "NH2", "O21", 2.58)),
        "reactive": 12.70,
    },
    {
        "key": "ethanol50", "panel": "C", "title": "50% ethanol", "time": "141.8 ns",
        "pdb": "ethanol50_141p8ns.pdb", "residues": (64, 253, 254, 276),
        "contacts": ((253, "ND2", "O24", 2.68), (276, "O", "N36", 3.06)),
        "reactive": 4.45,
    },
]

HB_COLORS = ("hb_purple", "hb_blue")


def one_atom(selection):
    count = cmd.count_atoms(selection)
    if count != 1:
        raise ValueError(f"Expected exactly one atom; got {count}: {selection}")
    return cmd.get_atom_coords(selection)


def atom_distance(left, right):
    return dist(one_atom(left), one_atom(right))


def configure():
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_fog", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("antialias", 2)
    cmd.set("ambient", 0.42)
    cmd.set("specular", 0.28)
    cmd.set("shininess", 25)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("stick_radius", 0.17)
    cmd.set("dash_gap", 0.25)
    cmd.set("dash_length", 0.25)
    cmd.set("dash_width", 2.8)
    cmd.set("label_distance_digits", 2)
    cmd.set_color("protein_gray", [0.73, 0.78, 0.81])
    cmd.set_color("protein_carbon", [0.47, 0.53, 0.57])
    cmd.set_color("ligand_carbon", [0.10, 0.62, 0.59])
    cmd.set_color("cys_orange", [0.93, 0.56, 0.12])
    cmd.set_color("hb_purple", [0.48, 0.33, 0.70])
    cmd.set_color("hb_blue", [0.21, 0.48, 0.78])


def prepare_system(info):
    name = info["key"]
    path = HERE / info["pdb"]
    if not path.is_file():
        raise FileNotFoundError(path)
    cmd.load(str(path), name)
    protein = f"({name} and not resn LIG)"
    ligand = f"({name} and resn LIG)"
    if cmd.count_atoms(protein) != 5166 or cmd.count_atoms(ligand) != 42:
        raise ValueError(f"Unexpected atom counts for {name}")

    cmd.hide("everything", name)
    cmd.dss(protein)
    cmd.show("cartoon", protein)
    cmd.color("protein_gray", protein)
    cmd.set("cartoon_transparency", 0.16, name)

    residue_selection = " or ".join(f"({name} and chain A and resi {n} and not resn LIG)"
                                    for n in info["residues"])
    cmd.select(f"{name}_site", residue_selection)
    cmd.show("sticks", f"{name}_site and not hydro")
    cmd.color("protein_carbon", f"{name}_site and elem C")
    cmd.color("cys_orange", f"{name} and resi 64 and elem S")

    cmd.show("sticks", f"{ligand} and not hydro")
    cmd.color("ligand_carbon", f"{ligand} and elem C")
    cmd.color("marine", f"{name} and elem N and not hydro")
    cmd.color("tv_red", f"{name} and elem O and not hydro")
    cmd.set("valence", 1, ligand)

    # Individual, validated heavy-atom contacts. No distance search is used.
    for index, (resi, protein_atom, ligand_atom, expected) in enumerate(info["contacts"]):
        left = f"({name} and chain A and resi {resi} and not resn LIG and name {protein_atom})"
        right = f"({name} and resn LIG and name {ligand_atom})"
        measured = atom_distance(left, right)
        if abs(measured - expected) > 0.10:
            raise ValueError(f"Contact mismatch in {name}: {left}, {right}: {measured:.2f} Å")
        object_name = f"{name}_hb{index + 1}"
        cmd.distance(object_name, left, right)
        cmd.color(HB_COLORS[index], object_name)
        cmd.hide("labels", object_name)
        print(f"{name}: {resi} {protein_atom}--ZQG {ligand_atom}: {measured:.2f} A")

    sg = f"({name} and chain A and resi 64 and name SG)"
    carbon = f"({name} and resn LIG and name C28)"
    measured = atom_distance(sg, carbon)
    if abs(measured - info["reactive"]) > 0.10:
        raise ValueError(f"SG--C28 distance mismatch in {name}: {measured:.2f} Å")
    if measured < 6.0:
        object_name = f"{name}_sg_c28"
        cmd.distance(object_name, sg, carbon)
        cmd.color("cys_orange", object_name)
        cmd.hide("labels", object_name)
    else:
        print(f"{name}: SG--C28 {measured:.2f} A, displaced; no contact line drawn")


def show_only(info, distances):
    cmd.disable("all")
    name = info["key"]
    cmd.enable(name)
    if distances:
        for suffix in ("hb1", "hb2", "sg_c28"):
            if suffix != "sg_c28" or info["reactive"] < 6.0:
                cmd.enable(f"{name}_{suffix}")


def render():
    configure()
    for info in SYSTEMS:
        prepare_system(info)

    # All three PDBs were prealigned on 331 C-alpha positions. Capture each
    # camera once from the complete aligned set and reuse it without reorienting.
    cmd.enable("all")
    cmd.orient("water and not resn LIG")
    cmd.zoom("water and not resn LIG", 4)
    overview_view = cmd.get_view()

    pocket_atoms = " or ".join(f"({s['key']} and resn LIG) or ({s['key']} and resi 64 and name SG)"
                               for s in SYSTEMS)
    cmd.orient(pocket_atoms)
    cmd.zoom(pocket_atoms, 6)
    pocket_view = cmd.get_view()

    for info in SYSTEMS:
        name = info["key"]
        show_only(info, distances=False)
        cmd.set_view(overview_view)
        cmd.png(str(OUT / f"{name}_overview.png"), width=2400, height=1800, dpi=300, ray=1)

        show_only(info, distances=True)
        cmd.set_view(pocket_view)
        cmd.png(str(OUT / f"{name}_contacts.png"), width=2400, height=1800, dpi=300, ray=1)
        cmd.save(str(OUT / f"{name}_editable.pse"))
        print("Rendered:", name)

    make_contact_panel()
    print("Outputs:", OUT)


def make_contact_panel():
    try:
        from PIL import Image, ImageDraw, ImageFont
    except ImportError:
        print("Pillow unavailable: individual PNGs are complete; panel montage skipped")
        return

    font_path = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
    try:
        title_font = ImageFont.truetype(font_path, 54)
        text_font = ImageFont.truetype(font_path, 34)
        note_font = ImageFont.truetype(font_path, 26)
    except OSError:
        title_font = text_font = note_font = ImageFont.load_default()

    panel = Image.new("RGB", (7200, 2140), "white")
    draw = ImageDraw.Draw(panel)
    for index, info in enumerate(SYSTEMS):
        x = index * 2400
        rendered = Image.open(OUT / f"{info['key']}_contacts.png").convert("RGB")
        panel.paste(rendered, (x, 150))
        draw.text((x + 95, 45), f"{info['panel']}   {info['title']}  |  {info['time']}",
                  font=title_font, fill="#183442")
        lines = []
        for resi, protein_atom, ligand_atom, measured in info["contacts"]:
            aa = {64: "Cys", 75: "Tyr", 208: "Arg", 238: "Arg", 253: "Asn", 276: "Asn"}[resi]
            lines.append(f"{aa}{resi} {protein_atom} ··· ZQG {ligand_atom}   {measured:.2f} Å")
        for n, (line, color) in enumerate(zip(lines, ("#7453ac", "#386cbc"))):
            y = 1955 + 47*n
            draw.ellipse((x+88, y+7, x+113, y+32), fill=color)
            draw.text((x+135, y), line, font=text_font, fill="#263b48")
        statement = (f"Cys64 SG ··· C28   {info['reactive']:.2f} Å"
                     + ("   (displaced)" if info["reactive"] > 6 else ""))
        draw.text((x+92, 2050), statement, font=note_font, fill="#a56d22")
    panel.save(OUT / "three_state_contact_panel.png", dpi=(300,300))


render()
