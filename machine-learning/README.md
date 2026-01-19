# Machine Learning & Big Data Analytics 

<p align="center">
  <img src="https://img.shields.io/badge/Tools-Python%20%7C%20TensorFlow%20%7C%20Keras-green?style=for-the-badge" alt="Tools">
  <img src="https://img.shields.io/badge/Topics-Deep%20Learning%20%7C%20CV%20%7C%20Scraping-blue?style=for-the-badge" alt="Topics">
  <img src="https://img.shields.io/badge/University-Sapienza-darkred?style=for-the-badge" alt="Sapienza">
  <img src="https://img.shields.io/badge/Year-2024%2F2025-orange?style=for-the-badge" alt="Year">
</p>

<p align="center">
  <b>Hands-on notebooks: denoising autoencoders, CNNs, transfer learning, and data collection pipelines</b>
</p>

---

## Overview

This repository contains lab homeworks and a final exam notebook developed for the **Machine Learning & Big Data Analytics Laboratory** (Sapienza University of Rome, A.Y. 2024/2025). 

The work spans:
- **Representation learning** (autoencoders)
- **Image classification** (from scratch vs transfer learning)
- **Big data acquisition/processing** via web scraping + tabular filtering

---

## Learning Focus

Key practical skills exercised throughout the lab:

| Skill | Description |
|-------|-------------|
| **Deep Learning** | Building and training models in TensorFlow/Keras (Model API, layers, training loops) |
| **ML Workflows** | End-to-end pipelines: preprocessing, augmentation, evaluation, inference |
| **Data Pipelines** | Collecting data from web pages and structuring it into analyzable tables |

---

## Notebooks & Activities examples

### 1) Denoising Autoencoder (MNIST)

Implemented a **convolutional autoencoder** that learns to reconstruct clean digits from noisy inputs.

```
Input (noisy) → Encoder → Latent Space → Decoder → Output (clean)
```

**Pipeline:**
- MNIST loading + normalization to `[0, 1]`
- Additive Gaussian noise with clipping
- Encoder/Decoder architecture with `Conv2D`, pooling, and upsampling
- Training: Adam optimizer + binary cross-entropy loss

**Visualization:**

| Noisy Input | Reconstructed Output |
|:-----------:|:--------------------:|
| ![noisy](https://via.placeholder.com/64x64?text=Noisy) | ![clean](https://via.placeholder.com/64x64?text=Clean) |

```python
# Noise addition
noise_factor = 0.5
x_train_noisy = x_train + noise_factor * np.random.normal(loc=0.0, scale=1.0, size=x_train.shape)
x_train_noisy = np.clip(x_train_noisy, 0., 1.)
```

---

### 2) Final Exam: Flower Classification (Custom CNN vs DenseNet121)

Image classification task on a small dataset with **4 classes** and **320 total images**.

#### Dataset Split
| Set | Purpose |
|-----|---------|
| Training | Model learning |
| Validation | Hyperparameter tuning & early stopping |
| Test | Final evaluation on `test.jpg` |

#### Models Compared

| Model | Type | Description |
|-------|------|-------------|
| **Custom CNN** | From Scratch | 4+ convolutional layers, trained end-to-end |
| **DenseNet121** | Transfer Learning | ImageNet pretrained, fine-tuned on flowers |

#### Pipeline

```python
# Data augmentation with ImageDataGenerator
datagen = ImageDataGenerator(
    rescale=1./255,
    rotation_range=20,
    width_shift_range=0.2,
    height_shift_range=0.2,
    horizontal_flip=True,
    validation_split=0.2
)
```

#### Evaluation Metrics
- ✅ Accuracy/Loss curves (training vs validation)
- ✅ Confusion matrix for class separability
- ✅ Inference comparison on `test.jpg`

---

### 3) Web Scraping (Big Data Acquisition)

Data collection exercises using Python scraping libraries.

**Tools:**
```python
import requests
from bs4 import BeautifulSoup
import pandas as pd
```

**Tasks:**
| Task | Description |
|------|-------------|
| Product Extraction | Title, description, price from e-commerce pages |
| Property Parsing | Extract attributes (e.g., "weightable", "available") |
| Pagination Handling | Iterate over search result pages |
| Data Structuring | Build `pandas.DataFrame` from scraped data |
| Filtering | Apply conditions (e.g., `price < threshold`) |

**Example Output:**

| Product | Price | Weightable |
|---------|-------|------------|
| Item A  | €29.99 | Yes |
| Item B  | €15.50 | No |
| Item C  | €42.00 | Yes |

---

### 4) ImageNet Inference (Pretrained Network)

Quick demo of applying a pretrained CNN for **top-k predictions**.

```python
from tensorflow.keras.applications import DenseNet121
from tensorflow.keras.applications.densenet import preprocess_input, decode_predictions

# Load and preprocess
img = load_img('test.jpg', target_size=(224, 224))
x = img_to_array(img)
x = preprocess_input(x)

# Predict
model = DenseNet121(weights='imagenet')
preds = model.predict(x)
print(decode_predictions(preds, top=5))
```

---

### 5) Cats vs Dogs (From Scratch)

Experiment outline for training an image classifier from scratch.

**Steps:**
1.  Download and unzip dataset
2.  Remove corrupted images (JFIF header check)
3.  Create datasets with `image_dataset_from_directory`
4.  Apply augmentation layers
5.  Train custom CNN architecture
6.  Evaluate on test set


---

## 🛠️ Tools & Stack

| Category | Tools |
|----------|-------|
| **Language** | Python  |
| **Deep Learning** | TensorFlow, Keras |
| **Data Processing** | NumPy, pandas |
| **Visualization** | Matplotlib, seaborn |
| **ML Utilities** | scikit-learn |
| **Web Scraping** | requests, BeautifulSoup4 |

---

## 📊 Key Results Summary

| Project | Metric | Result |
|---------|--------|--------|
| Denoising Autoencoder | Reconstruction Quality | Clean digit recovery from noisy input |
| Flower Classification (Custom) | Validation Accuracy | ~89% |
| Flower Classification (DenseNet) | Validation Accuracy | ~93% |
| Web Scraping | Records Collected | N products structured |

---
## 👤 Author

**Francesco Natali** (1945581)

</p>
