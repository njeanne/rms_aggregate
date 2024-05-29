#!/usr/bin/env python3

"""
Created on 18 Apr. 2023
"""

__author__ = "Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__email__ = "jeanne.n@chu-toulouse.fr"
__version__ = "1.2.0"

import argparse
import logging
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy
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


def extract_rmsd_data(conditions, method, out_dir, domain, md_time):
    """
    Extract the RMSD values of each sample and return the aggregated RMSD values for each frame.

    :param conditions: the conditions dataframe.
    :type conditions: pandas.DataFrame
    :param method: the method used for the aggregation.
    :type: str
    :param out_dir: the path to output directory.
    :type out_dir: str
    :param domain: the domain.
    :type domain: str
    :param md_time: the molecular dynamics simulation time.
    :type md_time: int
    :return: the aggregated data for each frame, the RMSD data grouped by condition and the conditions (in case one 
    condition is removed).
    :rtype: pandas.DataFrame, pandas.DataFrame, pandas.DataFrame
    """
    data_aggregated = {"frame": [], "condition": [], f"RMSD {method}": []}
    data_by_condition = {"sample": [], "condition": [], "RMSD": []}
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
        logging.info(f"Extracting data from {len(by_condition)} file{'s' if len(by_condition) > 1 else ''} for "
                     f"condition: {row_condition['condition']}")
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
            logging.info(f"{prefix:<80}{len(df_current['frames'])} frames.")
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

            # add the data by group for each sample
            data_by_condition["sample"] = data_by_condition["sample"] + [sample] * len(df_current["frames"])
            data_by_condition["condition"] = (data_by_condition["condition"] +
                                              [f"{row_condition['condition']} ({len(by_condition)})"] *
                                              len(df_current["frames"]))
            data_by_condition["RMSD"] = data_by_condition["RMSD"] + df_current["RMSD"].tolist()

        # aggregate
        data_aggregated["frame"] = data_aggregated["frame"] + frames
        data_aggregated["condition"] = data_aggregated["condition"] + [
            f"{row_condition['condition']} ({len(by_condition)})"] * len(frames)
        logging.info(f"\tComputing the RMSD {method} for "
                     f"{'the '+str(len(by_condition)) if len(by_condition) > 1 else 'the'} "
                     f"file{'s' if len(by_condition) > 1 else ''} of the condition: {row_condition['condition']}")
        aggregated = []
        for _, row_rmsd in df_raw.iterrows():
            if method == "median":
                aggregated.append(row_rmsd.median())
            elif method == "average":
                aggregated.append(row_rmsd.mean())
        data_aggregated[f"RMSD {method}"] = data_aggregated[f"RMSD {method}"] + aggregated
    # remove conditions if necessary
    for condition_to_remove in conditions_to_remove:
        conditions.drop(conditions[conditions["condition"] == condition_to_remove].index, inplace = True)

    df_by_condition = pd.DataFrame.from_dict(data_by_condition)
    out_by_condition = os.path.join(out_dir, f"RMSD_data_{domain.replace(' ', '-')}_{md_time}-ns.csv")
    df_by_condition.to_csv(out_by_condition, sep=',', index=False)
    logging.info(f"RMSD data written: {out_by_condition}")

    df_aggregated = pd.DataFrame.from_dict(data_aggregated)
    out_aggregated = os.path.join(out_dir, f"RMSD_aggregated_{method}_data_{domain.replace(' ', '-')}_{md_time}"
                                             f"-ns.csv")
    df_aggregated.to_csv(out_aggregated, sep=',', index=False)
    logging.info(f"RMSD {method} aggregated data written: {out_aggregated}")
    
    df_by_condition = pd.DataFrame.from_dict(data_by_condition)

    return df_aggregated, df_by_condition, conditions


def histogram_rmsd(src, src_ci, domain, md_time, dir_path, fmt, conditions_colors):
    """
    Create the histogram of the RMSD by condition.

    :param src: the RMSD data.
    :type src: pandas.DataFrame
    :param src_ci: the confidence intervals data.
    :type src_ci: pandas.DataFrame
    :param domain: the studied domain.
    :type domain: str
    :param md_time: the molecular dynamics duration.
    :type md_time: int
    :param dir_path: the output directory path.
    :type dir_path: str
    :param fmt: the plot output format.
    :type fmt: str
    :param conditions_colors: the colors of the conditions.
    :type conditions_colors: list

    """
    # clear the previous plot
    plt.clf()

    # create the histogram
    rmsd_ax = sns.histplot(data=src, x="RMSD", stat="density", kde=True, hue="condition", palette=conditions_colors,
                           alpha=0.5)
    # add the Means and Confidence intervals
    ci_y_coord = rmsd_ax.get_ylim()[1]
    for index, row in src_ci.iterrows():
        rmsd_ax.hlines(ci_y_coord, row["min confidence interval"], row["max confidence interval"], linestyle="--",
                       color=conditions_colors[index])
        rmsd_ax.plot(row["min confidence interval"], ci_y_coord, marker="|", color=conditions_colors[index])
        rmsd_ax.plot(row["mean"], ci_y_coord, marker="o", color=conditions_colors[index])
        rmsd_ax.plot(row["max confidence interval"], ci_y_coord, marker="|", color=conditions_colors[index])
        ci_y_coord = ci_y_coord + 0.05

    plot = rmsd_ax.get_figure()
    plt.suptitle(f"RMSD histogram on {md_time} ns: {domain}", fontsize="large", fontweight="bold")
    plt.xlabel(f"RMSD (\u212B)", fontweight="bold")
    plt.ylabel("Density", fontweight="bold")
    out_path_plot = os.path.join(dir_path, f"RMSD_histogram_{domain.replace(' ', '-')}_{md_time}-ns.{fmt}")
    plot.savefig(out_path_plot)
    logging.info(f"RMSD histogram by condition: {os.path.abspath(out_path_plot)}")



def data_for_stats(src):
    """
    Remove the number of samples in the conditions' column: "condition 1 (18)" -> "condition 1"

    :param src: the data.
    :type src: pandas.DataFrame
    :return: the
    """
    src_copy = src.copy(deep=True)
    src_copy["condition"] = src_copy["condition"].str.replace(r" \(\d+\)", "", regex=True)
    return src_copy


def write_stat_section_header(title):
    """
    Create the box for the statistics files.

    :param title: the title to write in the box.
    :type title: str
    :return: the box with the title.
    :rtype: str
    """

    box_size = len(title) + 2 * 6
    box = (f"{'#' * box_size}\n"
           f"#{' ' * (box_size - 2)}#\n"
           f"#{' ' * 5}{title}{' ' * 5}#\n"
           f"#{' ' * (box_size - 2)}#\n"
           f"{'#' * box_size}\n")
    return box


def data_description(data, padding, txt):
    """
    Create the description text.

    :param data: the data to describe.
    :type data: pandas.core.series.Series
    :param padding: the spaces between the key and the value of a result.
    :type padding: int
    :param txt: the text that presents the description.
    :type txt: str
    :return: the text of the data description.
    :rtype: str
    """
    description = scipy.stats.describe(data)
    txt = f"{txt}\n\t{'observations:':<{padding}}{description.nobs}"
    txt = f"{txt}\n\t{'minimum:':<{padding}}{description.minmax[0]}"
    txt = f"{txt}\n\t{'maximum:':<{padding}}{description.minmax[1]}"
    txt = f"{txt}\n\t{'mean:':<{padding}}{description.mean}"
    txt = f"{txt}\n\t{'variance:':<{padding}}{description.variance}"
    txt = f"{txt}\n\t{'skewness:':<{padding}}{description.skewness}"
    txt = f"{txt}\n\t{'kurtosis:':<{padding}}{description.kurtosis}\n"
    return txt


def compute_stats_by_condition(src, alpha_risk, padding_space, out_path):
    """
    Compute the statistics of the RMSD for all the samples by condition.

    :param src: the data for all the RMSD values.
    :type src: pandas.DataFrame
    :param alpha_risk: the alpha risk for the statistical tests.
    :type alpha_risk: float
    :param padding_space: the spaces between the key and the value of a result.
    :type padding_space: int
    :param out_path: the output path for the statistics results' file.
    :type out_path: str
    :return: the dataframe of the mean and confidence intervals by condition.
    :rtype: pandas.DataFrame
    """

    # get the conditions as a list, for loop performed to keep the conditions' order
    conditions = []
    for condition in src["condition"]:
        if condition not in conditions:
            conditions.append(condition)

    logging.info(f"Performing statistical tests for all the RMSD data for the following conditions "
                 f"\"{', '.join(conditions)}\":")
    txt_mean_ci = ""
    conf_intervals = {"condition": [], "mean": [], "min confidence interval": [], "max confidence interval": []}
    with open(out_path, "w") as out:
        title_box = write_stat_section_header("Data Description")
        out.write(title_box)
        for condition in conditions:
            condition_rows = src[f"RMSD"][src["condition"] == condition]
            data_description_text = data_description(condition_rows, padding_space,
                                                     f"\nDescription of the \"{condition}\" distribution:")
            out.write(data_description_text)
            txt_mean_ci = f"{txt_mean_ci}Condition \"{condition}\":\n"
            mean = numpy.mean(condition_rows)
            norm_ci = scipy.stats.norm.interval(1 - alpha_risk, loc=mean, scale=numpy.std(condition_rows))
            conf_intervals["condition"].append(condition)
            conf_intervals["mean"].append(mean)
            conf_intervals["min confidence interval"].append(norm_ci[0])
            conf_intervals["max confidence interval"].append(norm_ci[1])
            txt_mean_ci = f"{txt_mean_ci}\t{'Mean:':<{padding_space}}{mean}\n"
            txt_mean_ci = f"{txt_mean_ci}\t{'Confidence interval:':<{padding_space}}{norm_ci[0]} - {norm_ci[1]}\n\n"

        title_box = write_stat_section_header("Confidence Intervals")
        out.write(f"\n{title_box}\n")
        out.write(txt_mean_ci)

    logging.info(f"\tStatistical tests on all the samples by condition written: {out_path}")
    return pd.DataFrame.from_dict(conf_intervals)


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
    # clear the previous plot
    plt.clf()
    # create the lineplot
    rmsd_aggregated_ax = sns.lineplot(data=src, x="frame", y=f"RMSD {method}", hue="condition", palette=conditions_colors,
                           alpha=0.5)
    rmsd_aggregated_ax.set_xlim(min(src["frame"]), max(src["frame"]))
    plot = rmsd_aggregated_ax.get_figure()
    plt.suptitle(f"RMSD {method} on {md_time} ns: {domain}", fontsize="large", fontweight="bold")
    plt.xlabel("frame", fontweight="bold")
    plt.ylabel(f"RMSD {method} (\u212B)", fontweight="bold")
    out_path_plot = os.path.join(dir_path, f"RMSD_aggregated_{method}_{domain.replace(' ', '-')}_{md_time}-ns.{fmt}")
    plot.savefig(out_path_plot)
    logging.info(f"{method} RMSD aggregated by condition lineplot: {os.path.abspath(out_path_plot)}")


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
    # clear the previous plot
    plt.clf()
    # create the histogram
    rmsd_aggregated_ax = sns.histplot(data=src, x=f"RMSD {method}", stat="density", kde=True, hue="condition",
                           palette=conditions_colors, alpha=0.5)
    plot = rmsd_aggregated_ax.get_figure()
    plt.suptitle(f"RMSD {method} histogram on {md_time} ns: {domain}", fontsize="large", fontweight="bold")
    plt.xlabel(f"RMSD {method} (\u212B)", fontweight="bold")
    plt.ylabel("Density", fontweight="bold")
    out_path_plot = os.path.join(dir_path,
                                 f"RMSD_aggregated_{method}_histogram_{domain.replace(' ', '-')}_{md_time}-ns.{fmt}")
    plot.savefig(out_path_plot)
    logging.info(f"{method} RMSD aggregated by condition histogram: {os.path.abspath(out_path_plot)}")


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

def compute_stats_aggregated(src, method, alpha_risk, padding_space, out_path):
    """
    Compute the statistics of the aggregated data.

    :param src: the aggregated RMSD data.
    :type src: pandas.DataFrame
    :param method: the method used to gather the data.
    :type method: str
    :param alpha_risk: the alpha risk for the statistical tests.
    :type alpha_risk: float
    :param padding_space: the spaces between the key and the value of a result.
    :type padding_space: int
    :param out_path: the output path for the statistics results' file.
    :type out_path: str
    """

    # subsample by condition and add to a list for comparison tests
    aggregated_conditions_subsampling_list = []
    # get the conditions as a list, for loop performed to keep the conditions' order
    conditions = []
    for condition in src["condition"]:
        if condition not in conditions:
            conditions.append(condition)

    logging.info(f"Performing statistical tests on the aggregated data for the following conditions "
                 f"\"{', '.join(conditions)}\":")
    normality = True
    homoscedasticity = True
    data_description_text = ""
    with open(out_path, "w") as out:
        title_box = write_stat_section_header("Data Description")
        out.write(title_box)
        for condition in conditions:
            aggregated_condition_rows = src[f"RMSD {method}"][src["condition"] == condition]
            aggregated_conditions_subsampling_list.append(aggregated_condition_rows)

            if data_description_text == "":
                data_description_text = f"Description of the aggregated \"{condition}\" distribution:"
            data_description_text = data_description(aggregated_condition_rows, padding_space,
                                                     data_description_text)
            out.write(data_description_text)
            data_description_text = ""

            out.write(f"\nNormality test (Shapiro) for the \"{condition}\" distribution:\n")
            if len(aggregated_condition_rows) > 5000:
                logging.warning("\tShapiro Test (Normality) p-value may not be accurate as N > 5000")
                out.write("\t/!\ p-value may not be accurate as N > 5000")
            shapiro_test = scipy.stats.shapiro(aggregated_condition_rows)
            out.write(get_stat_results_as_text(shapiro_test, padding_space))
            message = (f"{'Interpretation:':<{padding_space}}\"{condition}\" distribution "
                       f"{'do not f' if shapiro_test.pvalue <= alpha_risk else 'f'}ollows a Normal law at \u03B1 risk "
                       f"of {alpha_risk} (p.value {'<=' if shapiro_test.pvalue <= alpha_risk else '>'} {alpha_risk})")
            out.write(f"\t{message}\n\n")
            if shapiro_test.pvalue < alpha_risk:
                normality = False

        out.write(f"\nVariance equality (Bartlett test) for the \"{', '.join(conditions)}\" distributions:\n")
        bartlett = scipy.stats.bartlett(*aggregated_conditions_subsampling_list)
        out.write(get_stat_results_as_text(bartlett, padding_space))
        message = (
            f"{'Interpretation:':<{padding_space}}For \"{', '.join(conditions)}\" distributions the variance are  "
            f"{'not' if bartlett.pvalue <= alpha_risk else ''} equals at \u03B1 risk "
            f"of {alpha_risk} (p.value {'<=' if bartlett.pvalue <= alpha_risk else '>'} {alpha_risk})")
        out.write(f"\t{message}\n\n")
        if bartlett.pvalue < alpha_risk:
            homoscedasticity = False

        title_box = write_stat_section_header("Statistical Tests")
        out.write(f"\n{title_box}")

        if normality and homoscedasticity:
            message = (
                " All the distributions follows a normal law and have equal variances, ANOVA and T-test will be "
                "applied. ")
            out.write(f"{message:*^150}\n")
            out.write(f"\nANOVA test for the \"{', '.join(conditions)}\" distributions:\n")
            test = scipy.stats.f_oneway(*aggregated_conditions_subsampling_list)
        else:
            if not normality and homoscedasticity:
                prefix_message = (
                    " At least one distribution do not follows a normal law and the distributions do not "
                    "have equal variances")
            elif not normality:
                prefix_message = " At least one distribution do not follows a normal law"
            else:
                prefix_message = " The distributions do not have equal variances"
            message = f"{prefix_message}, Kruskall-Wallis and Mann-Whitney U tests will be applied. "
            out.write(f"{message:*^150}\n")
            out.write(f"\nKruskall Wallis test for the \"{', '.join(conditions)}\" distributions:\n")
            test = scipy.stats.kruskal(*aggregated_conditions_subsampling_list)
        out.write(get_stat_results_as_text(test, padding_space))
        message = (
            f"{'Interpretation:':<{padding_space}}at least one distribution has a mean distinct from the others at a "
            f"\u03B1 risk of {alpha_risk} (p.value {'<=' if test.pvalue <= alpha_risk else '>'} {alpha_risk})")
        out.write(f"\t{message}\n\n")

        for i in range(len(conditions) - 1):
            distribution_i = src[f"RMSD {method}"][src["condition"] == conditions[i]]
            for j in range(i + 1, len(conditions)):
                distribution_j = src[f"RMSD {method}"][src["condition"] == conditions[j]]
                if normality:
                    out.write(
                        f"\nStudent test for the \"{conditions[i]}\" and \"{conditions[j]}\" distributions:\n")
                    test = scipy.stats.ttest_ind(distribution_i, distribution_j)
                else:
                    out.write(f"\nMann Withney U test for the \"{conditions[i]}\" and \"{conditions[j]}\" "
                              f"distributions:\n")
                    test = scipy.stats.mannwhitneyu(distribution_i, distribution_j)
                out.write(get_stat_results_as_text(test, padding_space))
                message = (f"{'Interpretation:':<{padding_space}}\"{conditions[i]}\" and \"{conditions[j]}\" have "
                           f"{'a different' if test.pvalue <= alpha_risk else 'the same'} mean at a \u03B1 risk of "
                           f"{alpha_risk} (p.value {'<=' if test.pvalue <= alpha_risk else '>'} {alpha_risk})")
                out.write(f"\t{message}\n\n")

    logging.info(f"\tStatistical tests on aggregated data written: {out_path}")


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
    parser.add_argument("-a", "--alpha", required=False, default=0.05, type=float,
                        help="the risk alpha to apply for the statistical tests.")
    parser.add_argument("-m", "--method-aggregation", required=False, default="median", choices=["median", "average"],
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
    parser.add_argument("-p", "--padding", required=False, default=40, type=int,
                        help="the space between the key and the value of a result in the statistical files.")
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
    rmsd_aggregated, rmsd_by_condition, data_conditions = extract_rmsd_data(data_conditions, args.method_aggregation,
                                                                            args.out, args.domain, args.md_time)

    # data for all the samples by condition
    rmsd_data_for_stats = data_for_stats(rmsd_by_condition)
    mean_conf_inter = compute_stats_by_condition(rmsd_data_for_stats, args.alpha, args.padding,
                                                 os.path.join(args.out,
                                                              f"statistics_RMSD_{args.domain.replace(' ', '-')}_"
                                                              f"{args.md_time}-ns.txt"))
    histogram_rmsd(rmsd_by_condition, mean_conf_inter, args.domain, args.md_time, args.out, args.format,
                   data_conditions["color"].to_list())


    # aggregated data by condition
    lineplot_aggregated_rmsd(rmsd_aggregated, args.md_time, args.out, args.format, data_conditions["color"].to_list(),
                             args.method_aggregation, args.domain)
    histogram_aggregated_rmsd(rmsd_aggregated, args.md_time, args.out, args.format,
                              data_conditions["color"].to_list(), args.method_aggregation, args.domain)
    rmsd_data_aggregated_for_stats = data_for_stats(rmsd_aggregated)
    compute_stats_aggregated(rmsd_data_aggregated_for_stats, args.method_aggregation, args.alpha, args.padding,
                             os.path.join(args.out, f"statistics_{args.method_aggregation}_aggregated_RMSD_"
                                                    f"{args.domain.replace(' ', '-')}_{args.md_time}-ns.txt"))