<div align="center">

# 🤖 Machine Learning & Big Data Analytics

**Tabular prediction · Computer vision · Web scraping**

<img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python">
<img src="https://img.shields.io/badge/TensorFlow-FF6F00?style=for-the-badge&logo=tensorflow&logoColor=white" alt="TensorFlow">
<img src="https://img.shields.io/badge/scikit--learn-F7931E?style=for-the-badge&logo=scikit-learn&logoColor=white" alt="scikit-learn">
<img src="https://img.shields.io/badge/Sapienza-8b1a1a?style=for-the-badge" alt="Sapienza">

<sub>MSc in Statistical Methods and Applications · Sapienza University of Rome</sub>

</div>

---

## 📋 What is here

Sixteen notebooks and scripts from the Machine Learning and Big Data Analytics labs, grouped into three tracks: predicting a binary outcome from tabular data, image classification with convolutional networks, and collecting data from the web.

One result runs through the vision work and is worth stating up front. A convolutional network trained from scratch on 320 images reached **28% validation accuracy** while hitting 84% on training data. The same task, handed to a DenseNet121 pretrained on ImageNet, reached **97%**. The gap between those two numbers is the entire argument for transfer learning on small datasets, and it shows up twice in this folder.

| Track | Files | Data included |
|---|---|---|
| 📊 **Tabular** | 4 notebooks | ✅ `adult_census_income.csv` (32,561 rows) |
| 🖼️ **Vision** | 7 notebooks | ❌ image archives not redistributed |
| 🕸️ **Scraping** | 5 files | live websites |

---

## 📊 Track 1 — Income prediction

**Task:** predict whether a person earns more than $50K a year from census demographics. The dataset is the UCI Adult census extract, 32,561 records and 14 predictors covering age, education, occupation, marital status, hours worked and capital gains.

The raw `occupation` column has 14 levels plus missing values coded as `?`. Both get handled before modelling: rows with unknown occupation are dropped, and the remaining categories collapse into five broader groups.

```python
occupation_map = {
    'Adm-clerical': 'Office & Tech Roles',  'Exec-managerial': 'Office & Tech Roles',
    'Prof-specialty': 'Office & Tech Roles', 'Tech-support': 'Office & Tech Roles',
    'Sales': 'Sales & Personal Services',   'Other-service': 'Sales & Personal Services',
    'Craft-repair': 'Technical & Mechanical Work',
    'Machine-op-inspct': 'Technical & Mechanical Work',
    'Protective-serv': 'Security & Protection',
    'Armed-Forces': 'Agriculture & Military',
}
df = df[df['occupation'] != '?']
```

### Results

Cross-validated scores after `GridSearchCV` tuning, 70/30 split:

| Model | Score | Notebook |
|---|---|---|
| Logistic regression | 0.846 | [`tabular_02`](tabular_02_income_prediction_classical_models.ipynb) |
| Random forest | 0.860 | [`tabular_02`](tabular_02_income_prediction_classical_models.ipynb) |
| **Gradient boosting** | **0.865** | [`tabular_02`](tabular_02_income_prediction_classical_models.ipynb) |
| Feed-forward network, best of 16 configs | 0.852 | [`tabular_03`](tabular_03_income_prediction_neural_network.ipynb) |

Held-out accuracy lands between 0.85 and 0.87 for every model. Gradient boosting takes it by about two points over logistic regression, and the best of sixteen tuned neural networks beats neither. On structured data of this size and shape, the extra capacity buys nothing.

<details>
<summary><b>The neural network, and why the search was worth running</b></summary>

<br>

[`tabular_03`](tabular_03_income_prediction_neural_network.ipynb) puts the numeric columns through `StandardScaler` and the categorical ones through `OneHotEncoder`, wired together in a `ColumnTransformer` so the fit happens on training data only. It then sweeps 16 configurations with `ParameterGrid`:

```python
param_grid = {
    'learning_rate': [0.001, 0.01],
    'dropout_rate' : [0.2, 0.3],
    'units'        : [[64, 32, 16, 8], [128, 64, 32, 16]],
    'batch_size'   : [32, 64]
}
early_stop = EarlyStopping(patience=5, restore_best_weights=True, monitor='val_loss')
```

Each configuration trains for up to 50 epochs with early stopping. All 16 land between **0.8460 and 0.8516** — a spread of six tenths of a percentage point across two learning rates, two dropout rates, two depths and two batch sizes.

That flatness is the finding. When the architecture stops mattering, the ceiling is in the data rather than the model, which is exactly what the tree ensembles had already suggested.

</details>

Also in this track: [`tabular_04`](tabular_04_binary_classification_mlp.ipynb), a small MLP (32→16→1) on an anonymous 30-feature dataset, evaluated by accuracy, sensitivity and AUC. It reaches an **AUC of 0.998**, which is high enough to be worth treating with suspicion rather than pride: on a genuinely hard problem this does not happen, so the features are either near-deterministic or the split leaks.

---

## 🖼️ Track 2 — Computer vision

### The flower classification exam

**Task:** four flower species, 80 images each, 320 in total. Build a CNN with at least four convolutional layers, diagnose it, then repeat with a pretrained DenseNet121 and compare.

| Class | Folder | Species |
|---|---|---|
| Yellow flower | `0` | Daffodil |
| White flower | `2` | Lily-of-the-valley |
| Purple flower | `4` | Crocus |
| Sunflower | `9` | Sunflower |

### Results

| Model | Training acc. | Validation acc. | Epochs |
|---|---|---|---|
| Custom CNN | 0.838 | **0.281** | 15 |
| Custom CNN + augmentation + BatchNorm | 0.800 | **0.563** | 4 |
| **DenseNet121** (ImageNet weights) | 0.936 | **0.969** | 10 |

The first model is a textbook overfit. Training accuracy climbs past 84% while validation accuracy wanders between 25% and 62% and finishes at 28%, barely above the 25% you would get by guessing among four classes. The confusion matrix puts the errors mostly between crocus and sunflower against lily-of-the-valley.

Heavier augmentation and batch normalisation roughly double validation accuracy to 56%, which is progress and still not usable.

DenseNet121 with frozen ImageNet features clears 95% within four epochs and settles at 96.9%. Eighty images per class is far too little to learn edge and texture detectors from nothing, and far more than enough to fine-tune detectors somebody else already learned.

On the held-out `test.jpg`, both models agree on **Purple Flower**: the custom CNN at 0.9721 confidence, DenseNet at 0.9923. Two architectures trained differently arriving at the same answer is reassuring, though a confident prediction from a model that is 28% accurate overall is a good reminder that confidence and correctness are different things.

### The same lesson, second dataset

[`vision_05`](vision_05_cnn_from_scratch_dogs_vs_carps.ipynb) trains a network from scratch to separate dogs from carps. Training accuracy 0.919, validation accuracy 0.666. Same shape, same cause.

<details>
<summary><b>Autoencoders and the rest of the vision track</b></summary>

<br>

| Notebook | What it does |
|---|---|
| [`vision_01`](vision_01_random_forest_tuning_mnist.ipynb) | Random forest on MNIST, sweeping `max_features` at fixed tree count |
| [`vision_02`](vision_02_denoising_autoencoder_mnist.ipynb) | Convolutional denoising autoencoder on MNIST |
| [`vision_03`](vision_03_denoising_autoencoder_cifar10.ipynb) | The same idea on CIFAR-10, so colour instead of greyscale |
| [`vision_04`](vision_04_reference_keras_cats_vs_dogs.ipynb) | The Keras reference pipeline that `vision_05` adapts, kept unmodified and credited |
| [`vision_06`](vision_06_pretrained_densenet_inference.ipynb) | Top-k inference with an off-the-shelf ImageNet DenseNet121 |

The denoising autoencoders add Gaussian noise to the input and train the network to reconstruct the clean image:

```python
noise_factor = 0.5
x_train_noisy = x_train + noise_factor * np.random.normal(0.0, 1.0, x_train.shape)
x_train_noisy = np.clip(x_train_noisy, 0., 1.)
```

Encoder and decoder are built with `Conv2D`, pooling and upsampling, trained with Adam and binary cross-entropy. Reconstruction quality is judged visually, side by side, since there is no accuracy metric to report for this task.

</details>

---

## 🕸️ Track 3 — Web scraping

Five exercises in pulling structured data out of pages that were not designed to give it up.

| File | Target | Technique |
|---|---|---|
| [`scraping_01`](scraping_01_product_page_parsing.py) | americanmadestore.us | Single product page: title, description, price, attributes |
| [`scraping_02`](scraping_02_search_and_pagination.py) | americanmadestore.us, national-hardware.com | Search results across multiple pages |
| [`scraping_03`](scraping_03_exercises.ipynb) | americanmadestore.us | Parsing exercises, results into a DataFrame |
| [`scraping_04`](scraping_04_selenium_subito.ipynb) | subito.it | Selenium, for a page that renders its listings in JavaScript |
| [`scraping_05`](scraping_05_exam_books_toscrape.ipynb) | books.toscrape.com | Graded exam: catalogue crawl with pagination |

```python
headers = {'User-Agent': 'frank/1.0'}
p = requests.get(url, headers=headers)
page = BeautifulSoup(p.content, 'lxml')
```

`requests` plus BeautifulSoup covers four of the five. Subito.it needs Selenium because the listings are not in the HTML that arrives from the server.

---

## ⚙️ Running the code

```
tensorflow  keras  scikit-learn  pandas  numpy  matplotlib  seaborn
requests  beautifulsoup4  lxml  selenium
```

> **Note on data.** The tabular notebooks read `adult_census_income.csv`, which is included, so those run as they are. The vision notebooks expect image archives (`jpg.zip`, `pets.zip`, `kagglecatsanddogs_5340.zip`) and `tabular_04` expects `dati30.csv`; these were provided by the course and are not redistributed here. Every notebook keeps its stored outputs, so the results and figures are readable without re-running anything.

Most of this was written in Colab. Paths have been made relative where the data ships with the repository.

---

## 🧰 Stack

| | |
|---|---|
| **Language** | Python |
| **Deep learning** | TensorFlow, Keras (Sequential and Model API) |
| **Classical ML** | scikit-learn: logistic regression, random forest, gradient boosting, `GridSearchCV` |
| **Data** | pandas, NumPy |
| **Plotting** | Matplotlib, seaborn |
| **Scraping** | requests, BeautifulSoup4, lxml, Selenium |
