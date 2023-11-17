import torch
import torchvision.transforms as transforms
from pathlib import Path
import numpy as np
from PIL import Image
import pytorch_msssim
import matplotlib.pyplot as plt

def load_image(img_path):
    img = Image.open(img_path)
    transform = transforms.ToTensor()
    return transform(img).unsqueeze(0)
    
def load_image_contrast(img_path):
    img = Image.open(img_path)
    transform = transforms.ToTensor()
    return transform(img).unsqueeze(0), np.array(img.convert('L'))  # Also return numpy array

def calculate_rms_contrast(image):
    mean_intensity = np.mean(image)
    rms_contrast = np.sqrt(np.mean((image - mean_intensity)**2))
    return rms_contrast

img_dir = Path('/home/alingold/pix2pix/results/combined_plant_0bp_noBlurryRemoval_oldtest/test_latest/images')
#img_dir = Path('/home/alingold/BCI/PyramidPix2pix/results/combined_plant_0bp_noBlurryRemoval_20231012_noblurry100/test_latest/images')
real_images = list(img_dir.glob('*real_B.png'))
fake_images = [img_dir / f"{img.stem.replace('_real_B', '_fake_B')}.png" for img in real_images]


mse_losses = []
ssim_values = []
cross_entropy_values = []
rms_contrasts = []

for real_img_path, fake_img_path in zip(real_images, fake_images):
    _, real_img_np = load_image_contrast(real_img_path)
    _, fake_img_np = load_image_contrast(fake_img_path)
    
    rms_contrast_real = calculate_rms_contrast(real_img_np)
    rms_contrast_fake = calculate_rms_contrast(fake_img_np)
    
    rms_contrasts.append((rms_contrast_real, rms_contrast_fake))

# Separate RMS contrasts for real and fake images
rms_contrasts_real, rms_contrasts_fake = zip(*rms_contrasts)

for real_img_path, fake_img_path in zip(real_images, fake_images):
    real_img = load_image(real_img_path)
    fake_img = load_image(fake_img_path)
    
    # Upscale the fake image to match the real image dimensions and convert to color
    #fake_img = upscale_image(fake_img, size=(256, 256))
    #fake_img = fake_img.repeat(1, 3, 1, 1)
    
    mse_loss = torch.nn.functional.mse_loss(real_img, fake_img)
    mse_losses.append(mse_loss.item())

    ssim_value_torch = pytorch_msssim.ssim(real_img, fake_img)  # PyTorch SSIM
    ssim_values.append(ssim_value_torch.item())

    cross_entropy = torch.nn.functional.binary_cross_entropy(real_img, fake_img)
    cross_entropy_values.append(cross_entropy.item())
    
    # ssim_values.append(ssim_value_torch)

average_mse = np.mean(mse_losses)
average_ssim = np.mean(ssim_values)
average_cross_entropy = np.mean(cross_entropy_values)

mae_losses = []

for real_img_path, fake_img_path in zip(real_images, fake_images):
    real_img = load_image(real_img_path)
    fake_img = load_image(fake_img_path)
    
    # Upscale the fake image to match the real image dimensions and convert to color
    #fake_img = upscale_image(fake_img, size=(256, 256))
    #fake_img = fake_img.repeat(1, 3, 1, 1)

    mae_loss = torch.nn.functional.l1_loss(real_img, fake_img)
    mae_losses.append(mae_loss.item())

average_mae = np.mean(mae_losses)

print(f"Average MAE: {average_mae}")
print(f"Average MSE: {average_mse}")
print(f"Average SSIM: {average_ssim}")
print(f"Average Cross Entropy: {average_cross_entropy}")
    
# Plotting
fig, axs = plt.subplots(2, 2, figsize=(15, 15))

# Assuming mse_losses, ssim_values, cross_entropy_values, and mae_losses are lists of equal length
axs[0, 0].scatter(rms_contrasts_real, mse_losses, label='Real Images')
axs[0, 0].set_xlabel('RMS Contrast')
axs[0, 0].set_ylabel('MSE Loss')
axs[0, 0].legend()

axs[0, 1].scatter(rms_contrasts_real, ssim_values, label='Real Images')
axs[0, 1].set_xlabel('RMS Contrast')
axs[0, 1].set_ylabel('SSIM Value')
axs[0, 1].legend()

axs[1, 0].scatter(rms_contrasts_real, cross_entropy_values, label='Real Images')
axs[1, 0].set_xlabel('RMS Contrast')
axs[1, 0].set_ylabel('Cross Entropy')
axs[1, 0].legend()

axs[1, 1].scatter(rms_contrasts_real, mae_losses, label='Real Images')
axs[1, 1].set_xlabel('RMS Contrast')
axs[1, 1].set_ylabel('MAE Loss')
axs[1, 1].legend()

plt.tight_layout()
plt.show()

