import os
import pandas as pd
import matplotlib.pyplot as plt
from src.helper_functions import Helper
helper = Helper()

root = "../../data/fatigue"
dir_data = root+"/plot_data"
dir_plots = root+"/plots"
rename_title = {
		"equivalent_load_range": "equivalent load range",
		"RootMEdg": "root edge-wise bending moment",
		"RootMFlp": "root flap-wise bending moment"
}
y_labels = {
		"equivalent_load_range": r"$\Delta M_{eq}$ (Nm)",
		"RootMEdg": r"partial damage $d$",
		"RootMFlp": r"partial damage $d$"
}

x_labels = {
		"equivalent_load_range": r"$N_{eq}$",
		"RootMEdg": r"U $(\frac{m}{s})$",
		"RootMFlp": r"U $(\frac{m}{s})$"
}
scale = {
		"equivalent_load_range": "log",
		"RootMEdg": "linear",
		"RootMFlp": "linear"
}
plt.rcParams.update({'font.size': 30})
for file in os.listdir(dir_data):
	plot_type = file[:file.find(".")]
	df = pd.read_csv(dir_data+"/"+file)
	columns = df.columns
	fig, ax = plt.subplots()
	for col in columns[1:]:
		ax.plot(df[columns[0]], df[col], label=col)
	helper.handle_axis(ax, legend=True, grid=True, x_label=x_labels[plot_type], y_label=y_labels[plot_type],
					   line_width=3, font_size=30, title=rename_title[plot_type], x_scale=scale[plot_type],
					   y_scale=scale[plot_type])
	if plot_type == "equivalent_load_range":
		handles, labels = ax.get_legend_handles_labels()
		order = [0, 2, 4, 1, 3, 5]  # reorder legend
		plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order])
	helper.handle_figure(fig, save_to=dir_plots+"/"+plot_type)
	
	
	
	
	
	
