from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import os
from tqdm import tqdm
from PIL import Image, ImageDraw, ImageFont


def mol_to_svg(mol_file_path, svg_file_path):
    mol = Chem.MolFromMolFile(mol_file_path)
    if not mol:
        print("Error: Couldn't read the molecule from the provided .mol file.")
        return

    # Using a molecular drawer
    drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)  # Here, 300x300 is the image size. You can adjust as needed.
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg_data = drawer.GetDrawingText().replace("svg:", "")

    with open(svg_file_path, "w") as svg_file:
        svg_file.write(svg_data)
        
def create_missing_image(svg_file_path):
    """Generate a gray SVG with the text 'Missing'."""
    width, height = 300, 300
    img = Image.new('RGBA', (width, height), (0, 0, 0, 0))  # Transparent background

    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype("arial.ttf", 40)  # Use a suitable font
    text_width, text_height = draw.textsize("Missing", font=font)
    draw.text(((width - text_width) / 2, (height - text_height) / 2), "Missing", fill="gray", font=font)

    with open(svg_file_path, "wb") as f:
        img.save(f, "SVG")

def process_compounds(base_path='comp_info'):
    for folder in tqdm(os.listdir(base_path)):
        folder_path = os.path.join(base_path, folder)
        if os.path.isdir(folder_path):
            mol_file = os.path.join(folder_path, f"{folder}.mol")
            svg_file = os.path.join(folder_path, f"{folder}_pic.svg")
            try:
                if os.path.exists(mol_file):
                    mol_to_svg(mol_file, svg_file)
                else:
                    create_missing_image(svg_file)
            except Exception as e:
                print(f"Error processing {folder}: {e}")
                create_missing_image(svg_file)

process_compounds()