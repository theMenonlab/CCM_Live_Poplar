import os
import re
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

# List of directories containing images to be averaged
#dirs = [
#    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_cannula_only_tile/test_latest/images",
#    "/home/alingold/BCI/PyramidPix2pix/results/combined_plant_0bp_noBlurryRemoval_cannula_only_tile/test_latest/images",
#    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_NoGan_cannula_only_tile/test_latest/images", # try with and without noGan
##    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_resnet_cannula_only/test_latest/images", # remove resnet from confidence map
#    "/home/alingold/pix2pix/results/plant_0bp_cannula_only_tile/test_latest/images" # blurry removal
#]

#dirs = [
#    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_cannula_only/test_latest/images",
#    "/home/alingold/BCI/PyramidPix2pix/results/combined_plant_0bp_noBlurryRemoval_cannula_only/test_latest/images",
#    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_NoGan_cannula_only/test_latest/images", # try with and without noGan
##    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_resnet_cannula_only/test_latest/images", # remove resnet from confidence map
#    "/home/alingold/pix2pix/results/plant_0bp_cannula_only/test_latest/images" # blurry removal
#]

dirs = [
    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_0929_x3/test_latest/images",
    "/home/alingold/BCI/PyramidPix2pix/results/combined_plant_0bp_noBlurryRemoval_0929_x3/test_latest/images",
    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_NoGan_0929_x3/test_latest/images", # try with and without noGan
#    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_resnet_cannula_only/test_latest/images", # remove resnet from confidence map
    "/home/alingold/pix2pix/results/plant_0bp_0929_x3/test_latest/images" # blurry removal
]

#dirs = [
#    "/home/alingold/pix2pix/results/plant_0bp_oldtest/test_latest/images",
#    "/home/alingold/BCI/PyramidPix2pix/results/combined_plant_0bp_noBlurryRemoval_oldtest/test_latest/images",
#    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_NoGan_oldtest/test_latest/images", # try with and without noGan
##    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_resnet_cannula_only/test_latest/images", # remove resnet from confidence map
#    "/home/alingold/pix2pix/results/plant_0bp_oldtest/test_latest/images"
#]

#dirs = [
#    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_20231012_noblurry100/test_latest/images",
#    "/home/alingold/BCI/PyramidPix2pix/results/combined_plant_0bp_noBlurryRemoval_20231012_noblurry100/test_latest/images",
#    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_NoGan_20231012_noblurry100/test_latest/images", # try with and without noGan
##    "/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_resnet_cannula_only/test_latest/images", # remove resnet from confidence map
#    "/home/alingold/pix2pix/results/plant_0bp_20231012_noblurry100/test_latest/images"
#]

# Output directory
output_dir = "/home/alingold/pix2pix/results/plant_0929_x3_log_confidence/test_latest/images"
os.makedirs(output_dir, exist_ok=True)

# Get list of image filenames from the first directory
image_filenames = os.listdir(dirs[0])

# List to store the average confidence for each image
avg_confidences = []

# Regex pattern for the required filenames
pattern = re.compile(r"\d+_fake_B\.png")

for image_filename in image_filenames:
    if not pattern.match(image_filename):
        continue
        
        
    # Use numpy to average the images
    imgs = [np.array(Image.open(os.path.join(dir_path, image_filename))) for dir_path in dirs]
    
    # Ensure that the image dimensions match across the directories
    assert all(img.shape == imgs[0].shape for img in imgs), f"Image shapes do not match for {image_filename}"
    
    # Average the images
    avg_img = np.mean(imgs, axis=0).astype(np.uint8)
    
    # Compute the confidence map and its average
    std_dev = np.std(imgs, axis=0)
    confidence_map = np.log(1 / (std_dev + 1))
    min_confidence = np.min(confidence_map)
    min_confidence = -3.9
    max_confidence = np.max(confidence_map)
    avg_confidence = np.mean(confidence_map)
    avg_confidences.append(avg_confidence)
    
    # Save the averaged image
    avg_img_pil = Image.fromarray(avg_img)
    avg_img_pil.save(os.path.join(output_dir, f"{image_filename}_conf_{avg_confidence:.3f}.png"))
    
    # Normalize the confidence map for display purposes
    confidence_map_norm = (confidence_map - min_confidence) / (max_confidence - min_confidence)
    fig, ax = plt.subplots()
    image = ax.imshow(confidence_map_norm, cmap='gray', vmin=0, vmax=1)  # Normalized to 0-1 for display

    
    # Remove the axes
    ax.axis('off')

    # Define the position of the colorbar
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    
    # Create a colorbar with specified ticks across the normalized range
    num_ticks = 5  # For example, 5 ticks
    tick_values = np.linspace(0, 1, num_ticks)
    tick_labels = [f"{v:.2f}" for v in np.linspace(min_confidence, max_confidence, num_ticks)]
    
    # Create the colorbar
    cbar = fig.colorbar(image, cax=cax, ticks=tick_values)
    
    # Manually set the height of the colorbar
    cbar.ax.set_aspect(18)  # The aspect ratio, increase or decrease for your desired height
    
    # Set the font size of the tick labels to be double the default size
    default_font_size = plt.rcParams.get('font.size')  # Get the default font size
    cbar.ax.set_yticklabels(tick_labels, fontsize=default_font_size * 2)  # Set text labels on the colorbar
    
    plt.savefig(os.path.join(output_dir, f"{image_filename}_confidence_map.png"), bbox_inches='tight', pad_inches=0)
    plt.close()

# Calculate and print the overall average confidence
overall_avg_confidence = np.mean(avg_confidences)
overall_std_confidence = np.std(avg_confidences)
print(f"Image averaging and confidence map generation completed.")
print(f"The overall average confidence for all images is: {overall_avg_confidence:.4f}")
print(f"The overall average confidence for all images is: {overall_std_confidence:.6f}")


