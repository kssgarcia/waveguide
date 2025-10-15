# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 16:35:26 2025

@author: agarz
"""

import cv2
import pytesseract
import pathlib
from pathlib import Path
from pdf2image import convert_from_path

pytesseract.pytesseract.tesseract_cmd = r"C:\Program Files\Tesseract-OCR\tesseract.exe"


# Read image
# easy_text_path = "easy_text.png"
easy_text_path = r"code_images\IMG_20250715_091209.jpg"
easy_img = cv2.imread(easy_text_path)
# Convert to text
text = pytesseract.image_to_string(easy_img)
print(text)

def image_to_text(input_path):
   """
   A function to read text from images.
   """
   img = cv2.imread(input_path)
   text = pytesseract.image_to_string(img)

   return text.strip()

# Define image path
# medium_text_path = "medium_text.png"

# # Extract text
# extracted_text = image_to_text(medium_text_path)
# print(extracted_text)

# def pdf_to_image(pdf_path, output_folder: str = "."):
#    """
#    A function to convert PDF files to images
#    """
#    # Create the output folder if it doesn't exist
#    if not Path(output_folder).exists():
#        Path(output_folder).mkdir()

#    pages = convert_from_path(pdf_path, output_folder=output_folder, fmt="png")

#    return pages

# pdf_path = "scanned_document.pdf"
# pdf_to_image(pdf_path, output_folder="documents")
