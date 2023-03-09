import pandas as pd
from matplotlib import colors
import matplotlib.pyplot as plt


def background_gradient(s: pd.Series, cmap="PuBu"):
    """generate a normalized color map for each column of the dataframe"""
    if s.min() > 0 and s.max() < 1:
        norm = colors.PowerNorm(2, vmin=0, vmax=1)
    else:
        norm = colors.PowerNorm(2, vmin=0, vmax=s.max())
    normed = norm(s.values)
    c = [colors.rgb2hex(x) for x in plt.colormaps.get_cmap(cmap)(normed)]
    return [f"background-color: {color}" for color in c]


def pretty_table(df: pd.DataFrame, cmap="PuBu"):
    return df.style.format(precision=3).apply(background_gradient, cmap=cmap)
