import pathlib
import os
import shutil
import matplotlib

class Helper():
    def __init__(self):
        pass

    def create_dir(self,
                   path_dir: str,
                   overwrite: bool = False,
                   add_missing_parent_dirs: bool = True,
                   raise_exception: bool = False,
                   print_message: bool = False) \
            -> tuple[str, bool]:
        return self._create_dir(path_dir, overwrite, add_missing_parent_dirs, raise_exception, print_message)

    @staticmethod
    def _create_dir(target: str,
                    overwrite: bool,
                    add_missing_parent_dirs: bool,
                    raise_exception: bool,
                    print_message: bool) -> tuple[str, bool]:
        msg, keep_going = str(), bool()
        try:
            if overwrite:
                if os.path.isdir(target):
                    shutil.rmtree(target)
                    msg = f"Existing directory {target} was overwritten."
                else:
                    msg = f"Could not overwrite {target} as it did not exist. Created it instead."
                keep_going = True
            else:
                msg, keep_going = f"Directory {target} created successfully.", True
            pathlib.Path(target).mkdir(parents=add_missing_parent_dirs, exist_ok=False)
        except Exception as exc:
            if exc.args[0] == 2:  # FileNotFoundError
                if raise_exception:
                    raise FileNotFoundError(f"Not all parent directories exist for directory {target}.")
                else:
                    msg, keep_going = f"Not all parent directories exist for directory {target}.", False
            elif exc.args[0] == 17:  # FileExistsError
                if raise_exception:
                    raise FileExistsError(f"Directory {target} already exists and was not changed.")
                else:
                    msg, keep_going = f"Directory {target} already exists and was not changed.", False
        if print_message:
            print(msg)
        return msg, keep_going

    @staticmethod
    def handle_figure(figure: matplotlib.figure,
                      file_figure: str=False,
                      show: bool=False,
                      size: tuple=(18.5, 10),
                      inches: int=100,
                      tight_layout: bool=True,
                      close: bool=True) -> matplotlib.figure:
        figure.set_size_inches(size)
        figure.set_dpi(inches)
        figure.set_tight_layout(tight_layout)
        if file_figure:
            figure.savefig(file_figure)
        if show:
            figure.show()
        if close: matplotlib.pyplot.close(figure)
        else: return figure

    @staticmethod
    def handle_axis(axis: matplotlib.pyplot.axis or list[matplotlib.pyplot.axis],
                    title: str=None or str or list,
                    grid: bool=False,
                    legend: bool=False,
                    legend_loc: int=0,
                    legend_columns: int=1,
                    x_label: str=False,
                    y_label: str or list[str]=False ,
                    z_label: str=False,
                    x_scale: str="linear",
                    y_scale: str="linear",
                    x_lim: tuple=None,
                    y_lim: tuple=None,
                    font_size: int=False,
                    line_width: int=None) -> matplotlib.pyplot.axis:
        axis = axis if type(axis) == list else [axis]
        y_label = y_label if type(y_label) == list else [y_label]
        for i, ax in enumerate(axis):
            if type(title) == str:
                if i == 0:
                    ax.set_title(title)
            else:
                if len(title) == 1:
                    ax.set_title(title[0])
                else:
                    ax.set_title(title[i])
            if x_lim is not None: ax.set_xlim(x_lim)
            if y_lim is not None: ax.set_ylim(y_lim)
            if font_size:
                ax.title.set_fontsize(font_size)
                ax.xaxis.label.set_fontsize(font_size)
                ax.yaxis.label.set_fontsize(font_size)
                if ax.name == "3d":
                    ax.zaxis.label.set_fontsize(font_size)
                ax.tick_params(axis='both', labelsize=font_size)
            if grid:
                axis[0].grid(grid)
            if x_label:
                ax.set_xlabel(x_label, labelpad=5)
            if y_label:
                ax.set_ylabel(y_label[i], labelpad=5)
            if z_label:
                ax.set_zlabel(z_label, labelpad=0)
            ax.set_xscale(x_scale)
            ax.set_yscale(y_scale)

        if legend or line_width:
            lines, labels = list(), list()
            for ax in axis:
                lines += ax.get_lines()
            if line_width:
                for line in lines:
                    line.set_linewidth(line_width)
            labels = [line.get_label() for line in lines if line.get_label()[0] != "_"]
            lines = [line for line in lines if line.get_label()[0] != "_"]
            if legend:
                if font_size:
                    axis[0].legend(lines, labels, ncol=legend_columns, prop={"size": font_size}, loc=legend_loc)
                else:
                    axis[0].legend(lines, labels, ncol=legend_columns, loc=legend_loc)
        return axis