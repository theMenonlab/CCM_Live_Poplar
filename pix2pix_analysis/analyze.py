import torch
import torchvision.transforms as transforms
from pathlib import Path
import numpy as np
from PIL import Image
import pytorch_msssim

def load_image(img_path):
    img = Image.open(img_path)
    transform = transforms.ToTensor()
    return transform(img).unsqueeze(0)
    
def upscale_image(image, size=(256, 256)):
    return torch.nn.functional.interpolate(image, size=size, mode='bilinear', align_corners=False)


img_dir = Path('/home/alingold/pix2pix/results/optrode_plant_ruipeng_nogan/test_latest/images')
#img_dir = Path('/home/alingold/BCI/PyramidPix2pix/results/combined_plant_0bp_noBlurryRemoval_20231012_noblurry100/test_latest/images')
real_images = list(img_dir.glob('*real_B.png'))
fake_images = [img_dir / f"{img.stem.replace('_real_B', '_fake_B')}.png" for img in real_images]


mse_losses = []
ssim_values = []
cross_entropy_values = []

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
    
    ssim_values.append(ssim_value_torch)
    
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
