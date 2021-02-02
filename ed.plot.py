import matplotlib as mpl
import itertools
import matplotlib.pyplot as plt
import pandas as pd
import collections
import argparse
mpl.use("Agg")


class Manhattan:
    def __init__(self,
                 infile,
                 X,
                 Y,
                 chr_list,
                 numbers,
                 x_label,
                 y_label,
                 name,
                 chr_name=None,
                 graph_width=14,
                 graph_height=7):
        self.infile = infile
        self.number = numbers
        self.chr_list = chr_list
        self.X = X
        self.Y = Y
        self.chr_name = chr_name
        self.x_label = x_label
        self.y_label = y_label
        self.name = name
        self.colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                       '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        self.df = self.get_data()
        self.figure = plt.figure(figsize=(graph_width, graph_height), frameon=True)
        # Creating the ax and modify it
        self.ax = self.figure.add_subplot(211)

        #set the tickt position
        self.ax.xaxis.set_ticks_position("bottom")
        self.ax.yaxis.set_ticks_position("left")

        #set the x-tict label
        self.ax.set_xlabel(self.x_label)
        self.ax.set_ylabel(self.y_label)
        # Hide the right and top spines
        self.ax.spines["top"].set_visible(False)
        self.ax.spines["right"].set_visible(False)
        self.ax.spines["bottom"].set_visible(False)

        # outward by 10 points
        for loc, spine in self.ax.spines.items():
            spine.set_position(('outward', 1))
        # Only draw spine between the y-ticks
        self.ax.spines['left'].set_bounds(0, 1)

    def get_chrs(self, chr_name='CHROM'):
        dict_chr_maker_numbers = collections.defaultdict(int)
        for chr_name in self.df[chr_name]:
            dict_chr_maker_numbers[chr_name] += 1
        sorted_dict_chr_maker_numbers = sorted(dict_chr_maker_numbers.items(), key=lambda x: x[1], reverse=True)
        chr_names = [i[0] for i in sorted_dict_chr_maker_numbers[0:self.number]]
        return sorted(chr_names)

    def get_chr_list(self, chr_list_file):
        chr_list = []
        with open(chr_list_file) as f:
            for line in f:
                chr_list.append(line.strip())
        return chr_list

    def plot_all_chr(self,
                     axis_text_size=12,
                     chr_text_size=12,
                     chr_name='CHROM',
                     dpi=600,
                     point_size=2.5):
        # chromosome_box_color = "#E5E5E5"
        # even_chromosome_color = "#1874CD"
        # odd_chromosome_color = "#4D4D4D"
        # plt.style.use('ggplot')
        chr_names = self.get_chr_list(self.chr_list)
        colors_iter = itertools.cycle(self.colors)
        start_pos = 0
        ticks = []
        i = 0
        max_y = []
        for name in chr_names:
            i += 1
            df_chr = self.df[self.df[chr_name] == name]
            color = next(colors_iter)
            y = df_chr[self.Y]
            max_y.append(max(y))
            self.ax.plot(start_pos + df_chr[self.X],
                         y,
                         ms=point_size,
                         marker='o',
                         ls="None",
                         mfc=color,
                         mec=color,
                         alpha=1)
            min_pos = start_pos
            max_pos = max(df_chr[self.X]) + start_pos
            # if i % 2 == 1:
            #     self.ax.axvspan(xmin=min_pos, xmax=max_pos, color="#E5E5E5")
            ticks.append((max_pos + min_pos) / 2)
            start_pos = max_pos

        # Putting the xticklabels
        self.ax.set_xticks(ticks)
        self.ax.set_xticklabels(chr_names)
        for tick in self.ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(axis_text_size)
        for tick in self.ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(chr_text_size)
            # tick.label.set_color('black')
            tick.label.set_rotation(90)
            # tick.label.set_fontsize(12)
        # Saving or plotting the figure
        mpl.rcParams['savefig.dpi'] = dpi
        mpl.rcParams['ps.papersize'] = "auto"
        mpl.rcParams['savefig.orientation'] = "landscape"
        self.ax.set_ylim(0, max(max_y))
        plt.savefig('{}.pdf'.format(self.name), dpi=dpi)

    def get_data(self):
        df = pd.read_csv(self.infile, sep='\t')
        df['Pos'] = (df['BIN_START'] + df['BIN_END'])/2
        return df

    def draw_one_chr(self, chr_name='chr'):
        df_chr = self.df[self.df[chr_name] == self.chr_name]
        return df_chr


def get_args():
    parser = argparse.ArgumentParser(description="plot pi value in the chr")
    parser.add_argument("-i", "--infile", help="infile")
    parser.add_argument("-X", "--pos", help="the name of the X value", default='Pos')
    parser.add_argument("-Y", "--value", help="the name of the Y value")
    parser.add_argument("-xl", "--x_label", help="the label name of the x axis")
    parser.add_argument("-yl", "--y_label", help="the label name of the y axis")
    parser.add_argument("-n", "--number", type=int, help="the number of chr to plot")
    parser.add_argument("-name", "--name_prefix", help="the name of the out to save prefix")
    parser.add_argument("-chr_list", "--chr_list", help="the name of the chr to plot")
    return parser.parse_args()


def main():
    arsg = get_args()
    infile = arsg.infile
    pos = arsg.pos
    value = arsg.value
    x_label = arsg.x_label
    y_label = arsg.y_label
    number = arsg.number
    name = arsg.name_prefix
    chr_list = arsg.chr_list
    plt_info = Manhattan(infile=infile,
                         X=pos,
                         Y=value,
                         chr_list=chr_list,
                         numbers=number,
                         x_label=x_label,
                         y_label=y_label,
                         name=name)
    plt_info.plot_all_chr()


main()