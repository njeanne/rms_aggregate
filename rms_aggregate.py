#!/usr/bin/env python3

"""
Created on 18 Apr. 2023
"""

__author__ = "Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__email__ = "jeanne.n@chu-toulouse.fr"
__version__ = "1.1.0"

import argparse
import logging
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scipy
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
    Extract the conditions, the paths and the colors.

    :param path: the path to the CSV file.
    :type path: str
    :return: the conditions.
    :rtype: pd.DataFrame
    """
    df = pd.read_csv(path, sep=",", header=None)
    df.columns = ["condition", "path", "color"]
    return df


def aggregate_rmsd(conditions, method, data_out_path):
    """
    Extract the RMSD values of each sample and return the aggregated RMSD values for each frame.

    :param conditions: the conditions dataframe.
    :type conditions: pandas.DataFrame
    :param method: the method used for the aggregation.
    :type: str
    :param data_out_path: the path to output data file.
    :type data_out_path: str
    :return: the aggregated data for each frame and the conditions (in case one condition is removed).
    :rtype: pandas.DataFrame, pandas.DataFrame
    """
    data = {"frames": [], "conditions": [], f"RMSD {method}": []}
    pattern = re.compile(".+?_(.+)\\.csv")
    frames = []
    conditions_to_remove = []
    for _, row_condition in conditions.iterrows():
        df_raw = pd.DataFrame()
        by_condition = [fn for fn in os.listdir(row_condition["path"]) if
                        fn.startswith("RMSD") and not "histogram" in fn and fn.endswith(".csv")]
        if len(by_condition) == 0:
            conditions_to_remove.append(row_condition["condition"])
            logging.warning(f"Condition {row_condition['condition']}: no RMSD files, this condition is skipped.")
            continue
        logging.info(f"Aggregating {len(by_condition)} file{'s' if len(by_condition) > 1 else ''} data for condition: "
                     f"{row_condition['condition']}")
        for item in sorted(by_condition):
            match = pattern.search(item)
            sample = None
            if match:
                sample = match.group(1)
            else:
                logging.error(f"\tNo match between the pattern '{pattern.pattern}' and the file name {item}.")
                sys.exit(1)
            df_current = pd.read_csv(os.path.join(row_condition["path"], item), sep=",")
            prefix = f"\t\t- {item}:"
            logging.info(f"{prefix:<60}{len(df_current['frames'])} frames.")
            if not frames:
                frames = df_current["frames"].to_list()
            elif len(frames) != len(df_current["frames"]):
                logging.error(f"\tNumber of frames for {item}, {len(df_current['frames'])} frames, differs from the "
                              f"previous number of frames, {len(frames)} frames.")
                sys.exit(1)
            elif frames != df_current["frames"].to_list():
                logging.error(f"\tThe frames for {item}, differs from the frames of the previous samples. Check the "
                              f"mask used for the RMS computation.")
                sys.exit(1)
            df_raw[sample] = df_current["RMSD"]
        data["frames"] = data["frames"] + frames
        data["conditions"] = data["conditions"] + [f"{row_condition['condition']} ({len(by_condition)})"] * len(frames)
        logging.info(f"Computing RMSD {method} for condition: {row_condition['condition']}")
        aggregated = []
        for _, row_rmsd in df_raw.iterrows():
            if method == "median":
                aggregated.append(row_rmsd.median())
            elif method == "average":
                aggregated.append(row_rmsd.mean())
        data[f"RMSD {method}"] = data[f"RMSD {method}"] + aggregated
    # remove conditions if necessary
    for condition_to_remove in conditions_to_remove:
        conditions.drop(conditions[conditions["condition"] == condition_to_remove].index, inplace = True)

    df = pd.DataFrame.from_dict(data)
    df.to_csv(data_out_path, sep=',', index=False)
    logging.info(f"Aggregated RMSD {method} data by condition written: {data_out_path}")

    return df, conditions


def lineplot_aggregated_rmsd(src, md_time, dir_path, fmt, conditions_colors, method, domain):
    """
    Create the lineplot of the aggregated RMSD.

    :param src: the data source.
    :type src: pandas.DataFrame
    :param md_time: the molecular dynamics duration.
    :type md_time: int
    :param dir_path: the output directory path.
    :type dir_path: str
    :param fmt: the plot output format.
    :type fmt: str
    :param conditions_colors: the colors of the conditions.
    :type conditions_colors: list
    :param method: the method used for the aggregation.
    :type: str
    :param domain: the studied domain.
    :type domain: str
    """
    rmsd_ax = sns.lineplot(data=src, x="frames", y=f"RMSD {method}", hue="conditions", palette=conditions_colors,
                           alpha=0.5)
    plot = rmsd_ax.get_figure()
    plt.suptitle(f"RMSD {method} on {md_time} ns: {domain}", fontsize="large", fontweight="bold")
    plt.xlabel("frames", fontweight="bold")
    plt.ylabel(f"RMSD {method} (\u212B)", fontweight="bold")
    out_path_plot = os.path.join(dir_path, f"RMSD_{method}_{domain.replace(' ', '-')}_{md_time}-ns.{fmt}")
    plot.savefig(out_path_plot)
    logging.info(f"RMSD {method} lineplot by condition: {os.path.abspath(out_path_plot)}")


def histogram_aggregated_rmsd(src, md_time, dir_path, fmt, conditions_colors, method, domain):
    """
    Create the histogram of the aggregated RMSD.

    :param src: the data source.
    :type src: pandas.DataFrame
    :param md_time: the molecular dynamics duration.
    :type md_time: int
    :param dir_path: the output directory path.
    :type dir_path: str
    :param fmt: the plot output format.
    :type fmt: str
    :param conditions_colors: the colors of the conditions.
    :type conditions_colors: list
    :param method: the method used for the aggregation.
    :type: str
    :param domain: the studied domain.
    :type domain: str
    """
    # clear the previous RMSD line plot
    plt.clf()
    # create the histogram
    rmsd_ax = sns.histplot(data=src, x=f"RMSD {method}", stat="density", kde=True, hue="conditions",
                           palette=conditions_colors, alpha=0.5)
    plot = rmsd_ax.get_figure()
    plt.suptitle(f"RMSD {method} histogram on {md_time} ns: {domain}", fontsize="large", fontweight="bold")
    plt.xlabel(f"RMSD {method} (\u212B)", fontweight="bold")
    plt.ylabel("Density", fontweight="bold")
    out_path_plot = os.path.join(dir_path, f"RMSD_{method}_histogram_{domain.replace(' ', '-')}_{md_time}-ns.{fmt}")
    plot.savefig(out_path_plot)
    logging.info(f"RMSD {method} histogram by condition: {os.path.abspath(out_path_plot)}")


def get_stat_results_as_text(stat_test, width):
    """
    Transform the test results as a string.

    :param stat_test: the statistical test.
    :type stat_test: scipy.stats
    :param width: the padding width.
    :type width: int
    :return: the string.
    :rtype: str
    """
    return f"\t{'p.value:':<{width}}{stat_test.pvalue}\n\t{'statistic:':<{width}}{stat_test.statistic}\n"


def compute_stats(src, method, out_path):
    """
    Compute the statistics of the distributions.

    :param src: the data.
    :type src: pandas.DataFrame
    :param method: the method used to gather the data.
    :type method: str
    :param out_path: the output path for the statistics results file.
    :type out_path: str
    """

    # remove the number of samples in the conditions' column: "condition 1 (18)" -> "condition 1"
    src["conditions"] = src["conditions"].str.replace(r" \(\d+\)", "", regex=True)

    # subsample by condition and add to a list for comparison tests
    conditions_subsampling_list = []
    conditions = list(set(src["conditions"].tolist()))

    logging.info(f"Performing statistical tests on the {method} for the following conditions "
                 f"\"{', '.join(conditions)}\":")
    padding = 20
    alpha_risk = 0.05
    with open(out_path, "w") as out:
        out.write("############################\n")
        out.write("#                          #\n")
        out.write("#     data description     #\n")
        out.write("#                          #\n")
        out.write("############################\n\n")
        normality = True
        variance_equality = True
        for condition in conditions:
            conditions_subsampling_list.append(src[f"RMSD {method}"][src["conditions"] == condition])
            condition_rows = src[f"RMSD {method}"][src["conditions"] == condition]
            out.write(f"Description of the \"{condition}\" distribution:\n")
            description = scipy.stats.describe(condition_rows)
            out.write(f"\t{'observations:':<{padding}}{description.nobs}\n")
            out.write(f"\t{'minimum:':<{padding}}{description.minmax[0]}\n")
            out.write(f"\t{'maximum:':<{padding}}{description.minmax[1]}\n")
            out.write(f"\t{'mean:':<{padding}}{description.mean}\n")
            out.write(f"\t{'variance:':<{padding}}{description.variance}\n")
            out.write(f"\t{'skewness:':<{padding}}{description.skewness}\n")
            out.write(f"\t{'kurtosis:':<{padding}}{description.kurtosis}\n")

            out.write(f"\nNormality test (Shapiro) for the \"{condition}\" distribution:\n")
            shapiro_test = scipy.stats.shapiro(condition_rows)
            out.write(get_stat_results_as_text(shapiro_test, padding))
            message = (f"{'Interpretation:':<{padding}}\"{condition}\" distribution "
                       f"{'do not f' if shapiro_test.pvalue <= alpha_risk else 'f'}ollows a Normal law at \u03B1 risk "
                       f"of {alpha_risk} (p.value {'<=' if shapiro_test.pvalue <= alpha_risk else '>'} {alpha_risk})")
            out.write(f"\t{message}\n\n")
            if shapiro_test.pvalue < alpha_risk:
                normality = False

        out.write(f"\nVariance equality (Bartlett test) for the \"{', '.join(conditions)}\" distributions:\n")
        bartlett = scipy.stats.bartlett(*conditions_subsampling_list)
        out.write(get_stat_results_as_text(bartlett, padding))
        message = (f"{'Interpretation:':<{padding}}For \"{', '.join(conditions)}\" distributions the variance are  "
                   f"{'not' if bartlett.pvalue <= alpha_risk else ''} equals at \u03B1 risk "
                   f"of {alpha_risk} (p.value {'<=' if bartlett.pvalue <= alpha_risk else '>'} {alpha_risk})")
        out.write(f"\t{message}\n\n")
        if bartlett.pvalue < alpha_risk:
            variance_equality = False

        out.write("\n############################\n")
        out.write("#                          #\n")
        out.write("#     statistical tests    #\n")
        out.write("#                          #\n")
        out.write("############################\n\n")

        if normality and variance_equality:
            message = (" All the distributions follows a normal law and have equal variances, ANOVA and T-test will be "
                       "applied. ")
            out.write(f"{message:*^150}\n")
            out.write(f"\nANOVA test for the \"{', '.join(conditions)}\" distributions:\n")
            test = scipy.stats.f_oneway(*conditions_subsampling_list)
        else:
            if not normality and variance_equality:
                prefix_message = (" At least one distribution do not follows a normal law and the distributions do not "
                                  "have equal variances")
            elif not normality:
                prefix_message = " At least one distribution do not follows a normal law"
            else:
                prefix_message = " The distributions do not have equal variances"
            message = f"{prefix_message}, Kruskall-Wallis and Mann-Whitney U tests will be applied. "
            out.write(f"{message:*^150}\n")
            out.write(f"\nKruskall Wallis test for the \"{', '.join(conditions)}\" distributions:\n")
            test = scipy.stats.kruskal(*conditions_subsampling_list)
        out.write(get_stat_results_as_text(test, padding))
        message = (f"{'Interpretation:':<{padding}}at least one distribution has a mean distinct from the others at a "
                   f"\u03B1 risk of {alpha_risk} (p.value {'<=' if test.pvalue <= alpha_risk else '>'} {alpha_risk})")
        out.write(f"\t{message}\n\n")

        for i in range(len(conditions) - 1):
            distribution_i = src[f"RMSD {method}"][src["conditions"] == conditions[i]]
            for j in range(i+1, len(conditions)):
                distribution_j = src[f"RMSD {method}"][src["conditions"] == conditions[j]]
                if normality:
                    out.write(f"\nStudent test for the \"{conditions[i]}\" and \"{conditions[j]}\" distributions:\n")
                    test = scipy.stats.ttest_ind(distribution_i, distribution_j)
                else:
                    out.write(f"\nMann Withney U test for the \"{conditions[i]}\" and \"{conditions[j]}\" "
                              f"distributions:\n")
                    test = scipy.stats.mannwhitneyu(distribution_i, distribution_j)
                out.write(get_stat_results_as_text(test, padding))
                message = (f"{'Interpretation:':<{padding}}\"{conditions[i]}\" and \"{conditions[j]}\" have "
                           f"{'a different' if test.pvalue <= alpha_risk else 'the same'} mean at a \u03B1 risk of "
                           f"{alpha_risk} (p.value {'<=' if test.pvalue <= alpha_risk else '>'} {alpha_risk})")
                out.write(f"\t{message}\n\n")

    logging.info(f"\tStatistical tests results written: {out_path}")


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

    The outputs are a line plot and an histogram with the aggregated RMSD values by condition, a CSV file with the 
    aggregated RMSD values used to create the plots and a statistics file.
    """
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--out", required=True, type=str, help="the path to the output directory.")
    parser.add_argument("-t", "--md-time", required=True, type=int,
                        help="the molecular dynamics duration in nanoseconds.")
    parser.add_argument("-d", "--domain", required=True, type=str,
                        help="Free text to specify the domain on which the RMS was computed.")
    parser.add_argument("-a", "--aggregation", required=False, default="median", choices=["median", "average"],
                        help="the aggregation method to apply.")
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
    rmsd_by_condition, data_conditions = aggregate_rmsd(data_conditions, args.aggregation,
                                                        os.path.join(args.out,
                                                                     f"RMSD_data_{args.domain.replace(' ', '-')}_"
                                                                     f"{args.md_time}-ns.csv"))
    lineplot_aggregated_rmsd(rmsd_by_condition, args.md_time, args.out, args.format, data_conditions["color"].to_list(),
                             args.aggregation, args.domain)
    histogram_aggregated_rmsd(rmsd_by_condition, args.md_time, args.out, args.format,
                              data_conditions["color"].to_list(), args.aggregation, args.domain)

    compute_stats(rmsd_by_condition, args.aggregation, os.path.join(args.out,
                                                                f"statistics_RMSD_data_{args.domain.replace(' ', '-')}"
                                                                f"_{args.md_time}-ns.txt"))