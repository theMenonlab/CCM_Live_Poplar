import numpy as np
import cv2
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
import numpy.ma as ma
import os
import scipy.io as sio
import skimage.transform as skt
from PIL import Image

def apply_circular_mask(image, center_x, center_y, mask_diameter):
    rows, cols = image.shape
    mask = np.zeros((rows, cols), dtype=np.uint8)
    cv2.circle(mask, (center_x, center_y), int(mask_diameter / 2), 1, -1)
    masked_image = ma.masked_array(image, mask=1 - mask)
    min_value_within_mask = masked_image.min()
    max_value_within_mask = masked_image.max()
    image = (image - min_value_within_mask) / (max_value_within_mask - min_value_within_mask)
    masked_image = image * mask

    return masked_image


def normalize_images(img, black_point, ccm):
    ref_imgs_normalized = np.zeros_like(ref_imgs)
    if ccm == 1: # check if it is ccm
        #img[ 94, 233] = 0.1 # fix broken pixel 92, 232#        img[ 93, 233] = 0.2 # fix broken pixel 92, 232
        #img[ 92, 232] = 0.1 # fix broken pixel 92, 232
        #img[ 93, 232] = 0.1 # fix broken pixel 92, 232
        #img[ 92, 231] = 0.1 # fix broken pixel 92, 232
        #img[ 93, 231] = 0.1 # fix broken pixel 92, 232
        #img[ 93, 230] = 0.1 # fix broken pixel 92, 232#        img[ 93, 233] = 0.2 # fix broken pixel 92, 232
        img[ 89, 219] = 0.1 # fix broken pixel 92, 232
        img[ 90, 220] = 0.1 # fix broken pixel 92, 232
        img[ 90, 219] = 0.1 # fix broken pixel 92, 232
        img[ 89, 220] = 0.1 # fix broken pixel 92, 232
    # Flatten the image to 1D array
    flat = img.flatten()
    # Sort the pixel intensities in ascending order
    sorted_pixels = np.sort(flat)
    # Compute the index of the black point
    black_index = int(black_point * len(sorted_pixels))
    # Set the black point as the intensity value below which all pixels will be considered black
    black_intensity = sorted_pixels[black_index]
    # Shift the intensity histogram to make the black point the new minimum
    img_shifted = img - black_intensity
    # Clip the intensities to the range [0, 1]
    img_shifted = np.clip(img_shifted, 0, None)
    # Normalize the intensities to the range [0, 1]
    img_min = img_shifted.min()
    img_max = img_shifted.max()
    ref_imgs_normalized = (img_shifted - img_min) / (img_max - img_min)
    return ref_imgs_normalized

def process_pix2pix(mnist_b, mnist_x):
        # Create the pix2pixData directory and the test subdirectory
    base_dir = r"C:\Users\menon\pytorch-CycleGAN-and-pix2pix\datasets\realtime"
    os.makedirs(base_dir, exist_ok=True)
    
    sub_dir = 'test'
    os.makedirs(os.path.join(base_dir, sub_dir), exist_ok=True)
    
    # Concatenate the images horizontally
    combined_image = np.concatenate((mnist_x, mnist_b), axis=1)

    # Scale the pixel values to the 16-bit range (0-65535)
    combined_image = (combined_image * 65535).astype(np.uint16)

    # Convert the combined image to a PIL Image object
    combined_image_pil = Image.fromarray(combined_image, mode='I;16')

    # Get the current date and time to include in the filename
    current_time = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

    # Save the combined image as n.png inside the test subdirectory
    # Add 1 to the index to start naming the files from 1.png
    file_path = os.path.join(base_dir, sub_dir, f'1_{current_time}.png')
    combined_image_pil.save(file_path)


reshape_shape = [256, 256]
# CCMcrop = [225, 480, 220, 475] 
CCMcrop = [237, 492, 223, 478] 
REFcrop = [430, 820, 425, 815]
center_x = 127
center_y = 128
mask_diameter = 240
black_point = 0

ref_imgs = []
ccm_imgs = []

Path = r"C:\Cannula Microscope\realtime\realtime.mat"
Struct = sio.loadmat(Path, verify_compressed_data_integrity=False)
Struct.keys()

ref_img = Struct["objRef"]
ccm_img = Struct["objCCM"]

ref_img = ref_img[REFcrop[2]:REFcrop[3], REFcrop[0]:REFcrop[1]]
ref_img = ref_img - np.min(ref_img)
ref_img = np.single(ref_img / np.max(ref_img))
ref_img = skt.resize(ref_img, reshape_shape)

ccm_img = ccm_img[CCMcrop[2]:CCMcrop[3], CCMcrop[0]:CCMcrop[1]]
ccm_img = ccm_img - np.min(ccm_img)
ccm_img = np.single(ccm_img / np.max(ccm_img))
ccm_img = skt.resize(ccm_img, reshape_shape)

ccm_img_pro = apply_circular_mask(ccm_img, center_x, center_y, mask_diameter)
ref_img_pro = normalize_images(ref_img, black_point, 0)
ccm_img_pro = normalize_images(ccm_img_pro, 0, 1)
print("Coordinates of max in ccm_img_pro:", np.unravel_index(np.argmax(ccm_img_pro), ccm_img_pro.shape))
print("max in ccm_img_pro:", np.max(ccm_img_pro))
process_pix2pix(ccm_img_pro, ref_img_pro)
