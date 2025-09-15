#!/usr/bin/env python3
import os
import argparse
import glob
import re
import tempfile
from reportlab.lib.pagesizes import letter, landscape
from reportlab.platypus import SimpleDocTemplate, Image, Paragraph, Table, TableStyle, PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from PIL import Image as PILImage


def create_comparison_slides(folder1, folder2, output_pdf):
    """
    Crea un PDF con diapositivas comparando imágenes de dos carpetas.
    
    Args:
        folder1: Ruta a la primera carpeta de imágenes
        folder2: Ruta a la segunda carpeta de imágenes
        output_pdf: Ruta del archivo PDF de salida
    """
    # Verificar que las carpetas existen
    if not os.path.isdir(folder1) or not os.path.isdir(folder2):
        raise ValueError("Ambas rutas deben ser directorios existentes")
    
    # Obtener lista de archivos de imagen en ambas carpetas
    ext_patterns = ["*.png", "*.jpg", "*.jpeg", "*.gif"]
    files1 = set()
    files2 = set()
    
    for pattern in ext_patterns:
        files1.update([os.path.basename(f) for f in glob.glob(os.path.join(folder1, pattern))])
        files2.update([os.path.basename(f) for f in glob.glob(os.path.join(folder2, pattern))])
    
    # Encontrar archivos comunes entre ambas carpetas
    common_files = sorted(list(files1.intersection(files2)))
    
    if not common_files:
        print(f"No se encontraron archivos de imagen comunes entre {folder1} y {folder2}")
        return
    
    # Crear documento PDF
    doc = SimpleDocTemplate(
        output_pdf,
        pagesize=landscape(letter),
        rightMargin=30,
        leftMargin=30,
        topMargin=30,
        bottomMargin=30
    )
    
    elements = []
    styles = getSampleStyleSheet()
    subtitle_style = styles["Heading2"]
    
    # Obtener nombres más descriptivos de carpetas (dos niveles)
    def get_descriptive_path(path):
        parts = os.path.normpath(path).split(os.sep)
        if len(parts) >= 2:
            return f"{parts[-2]}/{parts[-1]}"
        else:
            return parts[-1]
    
    folder1_name = get_descriptive_path(folder1)
    folder2_name = get_descriptive_path(folder2)
    
    # Función para redimensionar imágenes manteniendo la relación de aspecto
    def resize_image(img_path, max_width=380, max_height=280):
        img = PILImage.open(img_path)
        width, height = img.size
        
        # Calcular nueva dimensión manteniendo la relación de aspecto
        ratio = min(max_width/width, max_height/height)
        new_width = int(width * ratio)
        new_height = int(height * ratio)
        
        # Obtener la extensión del archivo original
        _, ext = os.path.splitext(img_path)
        
        # Crear un archivo temporal con la extensión correcta
        temp_file = tempfile.NamedTemporaryFile(suffix=ext, delete=False)
        temp_path = temp_file.name
        temp_file.close()
        
        # Redimensionar y guardar la imagen
        img = img.resize((new_width, new_height), PILImage.LANCZOS)
        img.save(temp_path)
        
        return temp_path, new_width, new_height
    
    # Para cada par de imágenes, crear una diapositiva en página completa
    for i, filename in enumerate(common_files):
        # Añadir salto de página excepto para la primera imagen
        if i > 0:
            elements.append(PageBreak())
        
        # Título de la imagen (eliminar extensión)
        image_name = os.path.splitext(filename)[0]
        # Mejorar el formato del título reemplazando _ por espacios y capitalizando palabras
        pretty_title = re.sub(r'[_-]', ' ', image_name).title()
        elements.append(Paragraph(f"<b>{pretty_title}</b>", subtitle_style))
        
        # Preparar imágenes
        img1_path = os.path.join(folder1, filename)
        img2_path = os.path.join(folder2, filename)
        
        # Asegurarse de que las imágenes existen
        if not os.path.exists(img1_path) or not os.path.exists(img2_path):
            print(f"Advertencia: No se puede encontrar {filename} en ambas carpetas. Saltando...")
            continue
        
        # Añadir timestamp único para evitar cacheo de imágenes
        img1_temp, w1, h1 = resize_image(img1_path)
        img2_temp, w2, h2 = resize_image(img2_path)
        
        # Crear objetos Image para ReportLab
        img1 = Image(img1_temp, width=w1, height=h1)
        img2 = Image(img2_temp, width=w2, height=h2)
        
        # Añadir leyendas para cada imagen
        data = [
            [folder1_name, folder2_name],
            [img1, img2]
        ]
        
        # Crear tabla con ancho completo de página
        t = Table(data, colWidths=[380, 380])
        t.setStyle(TableStyle([
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('BACKGROUND', (0, 0), (0, 0), colors.lightblue),
            ('BACKGROUND', (1, 0), (1, 0), colors.lightgreen),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
        ]))
        
        elements.append(t)
    
    # Generar el PDF
    doc.build(elements)
    print(f"PDF creado exitosamente: {output_pdf}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Crea un PDF con diapositivas comparando imágenes de dos carpetas")
    parser.add_argument("folder1", help="Ruta a la primera carpeta de imágenes")
    parser.add_argument("folder2", help="Ruta a la segunda carpeta de imágenes")
    parser.add_argument("--output", "-o", default="comparacion_imagenes.pdf", 
                        help="Ruta del archivo PDF de salida (por defecto: comparacion_imagenes.pdf)")
    
    args = parser.parse_args()
    
    try:
        create_comparison_slides(args.folder1, args.folder2, args.output)
    except Exception as e:
        print(f"Error: {e}")