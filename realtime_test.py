import os
import numpy as np
import time
import cv2
import subprocess
from options.test_options import TestOptions
from data import create_dataset
from models import create_model
from util.visualizer import save_images

# Define the path for the stopFlag
stopFlag_file = "stopFlag.txt"

print("Initializing...")

# Initialize test options
opt = TestOptions().parse()
opt.num_threads = 0
opt.batch_size = 1
opt.serial_batches = True
opt.no_flip = True
opt.display_id = -1

# Initialize dataset and model
print("Creating dataset and model...")
dataset = create_dataset(opt)
model = create_model(opt)
model.setup(opt)

# Create an output folder if it doesn't exist
output_folder = r"C:\Cannula Microscope\realtime_output"
file_name = "realtime_fake.png"
os.makedirs(output_folder, exist_ok=True)
print(f"Output will be saved to {output_folder}")

# Initialize processed files set
processed_files = set()

print("Starting monitoring...")


while True:
    # Check for the stopFlag
    if os.path.exists(stopFlag_file):
        print("Stop flag detected. Exiting...")
        break


    # Refresh dataset
    dataset = create_dataset(opt)

    for i, data in enumerate(dataset):
        if data['A_paths'][0] in processed_files:
            continue

        print(f"New file detected: {data['A_paths'][0]}")

        # Start the timer
        start_time = time.time()

        # Perform inference
        model.set_input(data)
        model.test()

        # Stop the timer and calculate the elapsed time
        elapsed_time = time.time() - start_time
        print(f"Inference time: {elapsed_time:.4f} seconds")

        # Get the base filename without extension
        base_filename = os.path.basename(data['A_paths'][0])
        base_filename_without_ext = os.path.splitext(base_filename)[0]

        # Create output filename
        output_filename = f"{base_filename_without_ext}_fake.png"
        
        # Save the result
        visuals = model.get_current_visuals()
        output_image_tensor = visuals['fake_B']
        output_image = (output_image_tensor.cpu().detach().numpy()[0, 0, :, :] + 1) / 2.0 * 255
        output_image = output_image.astype(np.uint8)
        output_image_path = os.path.join(output_folder, output_filename)
        cv2.imwrite(output_image_path, output_image)
        print(f"Saved to {output_image_path}")

        # Mark as processed
        processed_files.add(data['A_paths'][0])

    time.sleep(0.5)  # You can change this to be faster or slower