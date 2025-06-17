# SSCP 25 - Deep Learning for Cardiac Mechanics

This is a starter repository with 1000 sampled shape modes and blood pressures (Data/virtual_cohort_data.xlsx). The shape modes were sampled using a multivariate normal distribution sampler of the first 25 shape modes of the atlas. There are also two jupyter notebooks with some example ways of interacting with the sampled shape modes and the atlas.

## Set-up

1. Clone this github into your local device.

```
git clone https://github.com/bananna16/sscp25.git
```

2. Download the [EDES UKB atlas](https://drive.google.com/file/d/1nGlaEU_l6eJrSsk2DGX7Uqq0DgcDJU31/view?usp=sharing) into `sscp25/Data/`

3. Start-up the jupyter notebooks. I use [Visual Studio Code](https://code.visualstudio.com/) to interact with and use the jupyter notebooks. Of course, you can also [launch jupyter notebook through the terminal](https://docs.jupyter.org/en/latest/running.html) and intearct with them on the web browser.

4. Download the necessary imports. For example, `pip install numpy scipy pyvista`

- There's a special import that will be needed to interactively visualize the meshses uing pyvista: `pip install pyvista[jupyter]`

5. Run the jupyter notebooks:
- `VisualizeMeshes_AQ_2025_06_10.ipynb` has example code blocks and functions with how you can reconstruct the shapes from the z-scores and the atlas. There's also example code to visualize the shapes using pyvista.
- `GenUnloadedShape_AQ_2025_0518.ipynb` has code to generate volumes and masses from the shapes after they're reconstructed from the z-scores. We can also estimate the unloaded LV volume using an equation I generated after applying an adjusted [Garg et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35512290/) equation onto healthy individuals from UK Biobank. The equation in this notebook was just a linear equation fitted onto the data I had (with a r = 0.98). Then, I used scipy.optimize to adjust the PC scores until their volumes matched the target unloaded volume. There are two optimizes: one tries to keep the calculated LV mass the same where as the other one just optimizes on the LV volume. Trying to keep the LV mass constant is difficult when adjusting the meshes which causes the optimization to run for a while. These functions could very likely be improved upon a lot. There is also more code to visualize the unloaded shape.

If you have any questions, feel free to email me at anqi@ucsd.edu!
