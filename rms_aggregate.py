#!/usr/bin/env python3

"""
Created on 18 Apr. 2023
"""

__author__ = "Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__email__ = "jeanne.n@chu-toulouse.fr"
__version__ = "1.0.0"

import argparse
import logging
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def create_log(path, level):
    """Create the log as a text file and as a stream.

    :param path: the path of the log.
    :type path: str
    :param level: the level og the log.
    :type level: str
    :return: the logging:
    :rtype: logging
    """

    log_level_dict = {"DEBUG": logging.DEBUG,
                      "INFO": logging.INFO,
                      "WARNING": logging.WARNING,
                      "ERROR": logging.ERROR,
                      "CRITICAL": logging.CRITICAL}

    if level is None:
        log_level = log_level_dict["INFO"]
    else:
        log_level = log_level_dict[level]

    if os.path.exists(path):
        os.remove(path)

    logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s",
                        datefmt="%Y/%m/%d %H:%M:%S",
                        level=log_level,
                        handlers=[logging.FileHandler(path), logging.StreamHandler()])
    return logging


def get_conditions(path):
    """
    Extract the conditions, paths and colors.

    :param path: the path to the CSV file.
    :type path: str
    :return: the conditions.
    :rtype: pd.DataFrame
    """
    df = pd.read_csv(path, sep=",", header=None)
    df.columns = ["condition", "path", "color"]
    return df


def aggregate_rmsd(conditions):
    """
    Extract the RMSD values of each sample and return the median of the RMSD for each frame.

    :param conditions: the conditions dataframe.
    :type conditions: pd.DataFrame
    :return: the median for each frame.
    :rtype: pandas.DataFrame
    """
    data = {"frames": [], "conditions": [], "RMSD median": []}
    pattern = re.compile(".+?_(.+)\\.csv")
    frames = []
    for _, row_condition in conditions.iterrows():
        df_raw = pd.DataFrame()
        by_condition = [fn for fn in os.listdir(row_condition["path"]) if fn.startswith("RMSD") and fn.endswith(".csv")]
        logging.info(f"Aggregating {len(by_condition)} files data for condition: {row_condition['condition']}")
        for item in by_condition:
            logging.info(f"\t\t- {item}")
            match = pattern.search(item)
            if match:
                sample = match.group(1)
            else:
                logging.error(f"\tNo match between the pattern '{pattern.pattern}' and the file name {item}.")
                sys.exit(1)
            df_current = pd.read_csv(os.path.join(row_condition["path"], item), sep=",")
            if not frames:
                frames = df_current["frames"].to_list()
            elif len(frames) != len(df_current["frames"]):
                logging.error(f"\tNumber of frames for {item}, {len(df_current['frames'])} frames, differs from the "
                              f"previous number of frames, {len(frames)} frames.")
            elif frames != df_current["frames"].to_list():
                logging.error(f"\tThe frames for {item}, differs from the frames of the previous samples. Check the "
                              f"mask used for the RMS computation.")
            df_raw[sample] = df_current["RMSD"]
        data["frames"] = data["frames"] + frames
        data["conditions"] = data["conditions"] + [f"{row_condition['condition']} ({len(by_condition)})"] * len(frames)
        logging.info(f"Computing RMSD medians for condition: {row_condition['condition']}")
        medians = []
        for _, row_rmsd in df_raw.iterrows():
            medians.append(row_rmsd.median())
        data["RMSD median"] = data["RMSD median"] + medians
    return pd.DataFrame.from_dict(data)


def plot_aggregated_rmsd(src, md_time, dir_path, fmt, conditions_colors, region):
    """
    Plot the aggregated RMSD medians.

    :param src: the data source.
    :type src: pd.DataFrame
    :param md_time: the molecular dynamics duration.
    :type md_time: int
    :param dir_path: the output directory path.
    :type dir_path: str
    :param fmt: the plot output format.
    :type fmt: str
    :param conditions_colors: the colors of the conditions.
    :type conditions_colors: list
    :param region: the region.
    :type region: str
    """
    rmsd_ax = sns.lineplot(data=src, x="frames", y="RMSD median", hue="conditions", palette=conditions_colors,
                           alpha=0.5)
    plot = rmsd_ax.get_figure()
    plt.suptitle(f"RMSD: {region} {md_time} ns", fontsize="large", fontweight="bold")
    plt.xlabel("frames", fontweight="bold")
    plt.ylabel("RMSD medians (\u212B)", fontweight="bold")
    out_path_plot = os.path.join(dir_path, f"RMSD_{region.replace(' ', '-')}_{md_time}-ns.{fmt}")
    plot.savefig(out_path_plot)
    logging.info(f"RMSD median plot by condition: {os.path.abspath(out_path_plot)}")


if __name__ == "__main__":
    descr = f"""
    {os.path.basename(__file__)} v. {__version__}

    Created by {__author__}.
    Contact: {__email__}
    {__copyright__}

    Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or implied.

    Aggregate the RMS data in one plot to compare between various conditions. The RMS files are the RMSD CSV files from 
    the RMS analysis (https://github.com/njeanne/rms). The Molecular Dynamics simulation time must be identical.

    The input is a comma separated file without header which first column is the condition, the second column the path 
    of the directory containing the RMS analysis files and the third column the color in hexadecimal format. i.e:

    insertions,tests/inputs/insertions,#fc030b
    WT,tests/inputs/WT,#0303fc

    The output is a plot with the aggregated RMSD values by condition.
    """
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--out", required=True, type=str, help="the path to the output directory.")
    parser.add_argument("-t", "--md-time", required=True, type=int,
                        help="the molecular dynamics duration in nanoseconds.")
    parser.add_argument("-d", "--domain", required=True, type=str,
                        help="Free text to specify the domain on which the RMS was computed.")
    parser.add_argument("-x", "--format", required=False, default="svg",
                        choices=["eps", "jpg", "jpeg", "pdf", "pgf", "png", "ps", "raw", "svg", "svgz", "tif", "tiff"],
                        help="the output plots format: 'eps': 'Encapsulated Postscript', "
                             "'jpg': 'Joint Photographic Experts Group', 'jpeg': 'Joint Photographic Experts Group', "
                             "'pdf': 'Portable Document Format', 'pgf': 'PGF code for LaTeX', "
                             "'png': 'Portable Network Graphics', 'ps': 'Postscript', 'raw': 'Raw RGBA bitmap', "
                             "'rgba': 'Raw RGBA bitmap', 'svg': 'Scalable Vector Graphics', "
                             "'svgz': 'Scalable Vector Graphics', 'tif': 'Tagged Image File Format', "
                             "'tiff': 'Tagged Image File Format'. Default is 'svg'.")
    parser.add_argument("-l", "--log", required=False, type=str,
                        help="the path for the log file. If this option is skipped, the log file is created in the "
                             "output directory.")
    parser.add_argument("--log-level", required=False, type=str,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="set the log level. If the option is skipped, log level is INFO.")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("input", type=str,
                        help="the path to the CSV (comma separated without header) file which first column is the "
                             "condition, the second column the path of the directory containing the RMS analysis "
                             "files and the third column the color.")
    args = parser.parse_args()

    # create output directory if necessary
    os.makedirs(args.out, exist_ok=True)
    # create the logger
    if args.log:
        log_path = args.log
    else:
        log_path = os.path.join(args.out, f"{os.path.splitext(os.path.basename(__file__))[0]}.log")
    create_log(log_path, args.log_level)

    logging.info(f"version: {__version__}")
    logging.info(f"CMD: {' '.join(sys.argv)}")
    logging.info(f"MD simulation time: {args.md_time} ns")
    logging.info(f"Domain: {args.domain:>16}")

    data_conditions = get_conditions(args.input)
    rmsd_by_condition = aggregate_rmsd(data_conditions)
    plot_aggregated_rmsd(rmsd_by_condition, args.md_time, args.out, args.format, data_conditions["color"].to_list(),
                         args.domain)
