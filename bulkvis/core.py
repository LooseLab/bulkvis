"""core.py
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import traceback


def concat_files_to_df(file_list, **kwargs):
    """Return a pandas.DataFrame from a list of files
    Parameters
    ----------
    file_list : list
        List of files to be concatenated
    kwargs
        Any parameter used by pandas.read_csv except 'filepath_or_buffer'. These will be applied to all
        files in 'file_list'
    Returns
    -------
    pandas.DataFrame
    Raises
    ------
    pandas.errors.ParserError
        Raises pandas.errors.ParserError if input file(s) do not match expected format or shape.
    """
    kwargs = remove_kwargs(["filepath_or_buffer"], **kwargs)
    df_list = []
    for f in file_list:
        try:
            df_list.append(pd.read_csv(filepath_or_buffer=f, **kwargs))
        except pd.errors.ParserError as e:
            sys.exit(
                "ParserError\nUsually caused by an input file not being the expected format"
            )
        except Exception as e:
            traceback.print_exc()
            sys.exit(1)
    return pd.concat(df_list, ignore_index=True)


def remove_kwargs(remove_list, **kwargs):
    """Remove items from kwargs dict that may cause conflict with successive function calls"""
    # return {k: v for k, v in kwargs.items() if k not in remove_list}  # This iterates the entire dict
    for item in remove_list:  # This just iterates the remove_list
        _ = kwargs.pop(item, None)
    return kwargs


def length_stats(lengths):
    """Return count [COUNT], minimum [MIN], maximum [MAX], mean [MEAN], and N50 [N50] of an array
    Parameters
    ----------
    lengths : array_like
        List of integers
    Returns
    -------
    dictionary
    """
    return {
        "COUNT": int(len(lengths)),
        "MIN": int(np.min(lengths)),
        "MAX": int(np.max(lengths)),
        "MEAN": int(np.mean(lengths)),
        "N50": _get_n50(np.sort(lengths)),
    }


def _get_n50(lengths):
    """Return N50 statistic for a list of read lengths
    Parameters
    ----------
    lengths array_like
        List of sorted, ascending, integers
    Returns
    -------
    integer
    """
    return int(lengths[np.where(np.cumsum(lengths) >= np.sum(lengths) / 2)][0])


def readable_yield(num, suffix="B"):
    """Return a human readable string of yield using si/metric prefixes
    Parameters
    ----------
    num : int (or float)
        Integer of total number of bases
    suffix : str
        String to append to si/metric prefixes ['', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']
    Returns
    -------
    string
    """
    for unit in ["", "k", "M", "G", "T", "P", "E", "Z"]:
        if abs(num) < 1000:
            return "{n:3.2f} {u}{s}".format(n=num, u=unit, s=suffix)
        num /= 1000
    return "{n:3.1f} {u}{s}".format(n=num, u="Y", s=suffix)


def human_readable_yield(num: int, factor: int = 1000, suffix: str = "B") -> str:
    """Return a human readable string of a large number using SI unit prefixes
    Parameters
    ----------
    num : int
        A number to convert to decimal form
    factor : int
        The SI factor, use 1000 for SI units and 1024 for binary multiples
    suffix : str
        The suffix to place after the SI prefix, for example use B for SI units and iB for binary multiples
    Returns
    -------
    str
        Returns the input number formatted to two decimal places with the SI unit and suffix
    """
    for unit in ["", "k", "M", "G", "T", "P", "E", "Z"]:
        if abs(num) < factor:
            return "{n:3.2f} {u}{s}".format(n=num, u=unit, s=suffix)
        num /= factor
    return "{n:3.2f} {u}{s}".format(n=num, u="Y", s=suffix)


def top_n(df, field, n):
    """Print top N reads, by length, from a dataset
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing a Series with length values
    field : str
        The key for the Series containing length values
    n : int
        The number of values to print
    Returns
    -------
    None
    """
    df = df.sort_values(by=field, ascending=False)
    rows = df.filter([field], axis=1).head(n=n).reset_index()
    max_len = max(
        [len(str(r[field])) + len(str(r[field])) // 3 for i, r in rows.iterrows()]
    )
    for idx, row in rows.iterrows():
        print("{i}:\t{len:>{m},}".format(i=idx + 1, m=max_len, len=row[field]))
    return


def find_files_of_type(file_or_directory, file_extensions):
    """Return a list of pathlib.Path of files with chosen extensions
    Parameters
    ----------
    file_or_directory : str
        filepath or a directory
    file_extensions : list
        A list of lowercase file extensions including '.' e.g. ['.txt', '.paf']
    Returns
    -------
    list
        If files with extension are found return list of pathlib.Path objects, otherwise return empty list
    """
    file_or_directory = Path(file_or_directory).expanduser()
    if (
        file_or_directory.is_file()
        and "".join(file_or_directory.suffixes).lower() in file_extensions
    ):
        return [file_or_directory]
    elif file_or_directory.is_dir():
        return [
            x
            for x in file_or_directory.iterdir()
            if "".join(x.suffixes).lower() in file_extensions
        ]
    else:
        return []


def fuse_reads(seq_sum_df, paf_df, distance=10000, alt=True):
    """Find fused reads from sequencing_summary.txt and paf files
    Parse sequencing_summary.txt and mapping.paf files to infer reads that may
    have been incorrectly split my MinKNOW. This approach is based on read number
    and mapping, therefore a _good_ mapping is required.
    Parameters
    ----------
    seq_sum_df : pandas.DataFrame
        A pandas.DataFrame from a sequencing_summary.txt file, this must contain
        the columns `['channel', 'start_time', 'duration', 'run_id', 'read_id',
        'sequence_length_template', 'filename']`.
    paf_df : pandas.DataFrame
        A pandas.DataFrame from a .paf file, these are generated by minimap2.
        As this file type doesn't have headers the following parameters are
        the minimum required for using a file with this function
        `usecols=[0, 4, 5, 7, 8], names=['Qname', 'Strand', 'Tname', 'Tstart', 'Tend']`
    distance : int
        The distance, in bases, between the end coordinate of a read mapping and
        the start coordinate of successive read from the same channel. Defaults to 10000
    alt : bool
        Include alternate assemblies, default is True. If set to True (include
        alternate assemblies) the 'new' dataset may have more bases than the 'original'
        input dataset due to reads mapping to alternate contigs.
    Returns
    -------
    fused_reads_df : pandas.DataFrame
        pandas.DataFrame containing fused reads
    un_fused_reads_df : pandas.DataFrame
        pandas.DataFrame containing un-fused reads
    to_be_fused_reads_df : pandas.DataFrame
        pandas.DataFrame containing reads that are fused 'parts' in the same
        format as un_fused_reads_df
    """
    # TODO: raise error if required columns are not present
    # TODO: raise error if columns are not of correct types
    # Remove zero length reads and sort seq_sum_df
    seq_sum_df = seq_sum_df[seq_sum_df["sequence_length_template"] != 0].sort_values(
        by=["channel", "run_id", "start_time"]
    )
    # Create extra Series for finding fused reads
    seq_sum_df["next_read_id"] = seq_sum_df["read_id"].shift(-1)
    seq_sum_df["next_start_time"] = seq_sum_df["start_time"].shift(-1)
    seq_sum_df["next_end"] = seq_sum_df["next_start_time"] + seq_sum_df[
        "duration"
    ].shift(-1)
    seq_sum_df["next_sequence_length_template"] = seq_sum_df[
        "sequence_length_template"
    ].shift(-1)
    seq_sum_df["combined_length"] = (
        seq_sum_df["sequence_length_template"]
        + seq_sum_df["next_sequence_length_template"]
    )
    seq_sum_df["next_sequence_length_template"] = (
        seq_sum_df["next_sequence_length_template"].fillna(0).astype("int64")
    )
    seq_sum_df["combined_length"] = (
        seq_sum_df["combined_length"].fillna(0).astype("int64")
    )

    # Merge seq_sum_df and paf_df on read_id/Qname; this aligns the read and the mapping information
    df = pd.merge(seq_sum_df, paf_df, left_on="read_id", right_on="Qname", how="outer")
    # Merge df with paf_df on next_read_id/Qname; this aligns each read with it's
    # successor giving suffix '_A' and '_B' respectively
    df2 = pd.merge(
        df,
        paf_df,
        left_on="next_read_id",
        right_on="Qname",
        how="outer",
        suffixes=("_A", "_B"),
    )
    df2 = df2.dropna().reset_index()

    # If df2 had no rows, no merging has taken place
    if len(df2) == 0:
        return None, None, None

    # Condition where Qname (read_id) does NOT match
    not_qname = df2["Qname_A"] != df2["Qname_B"]
    # Condition where Strand matches
    yes_strand = df2["Strand_A"] == df2["Strand_B"]
    # Condition where Target (chromosome) matches
    yes_tname = df2["Tname_A"] == df2["Tname_B"]

    df2 = df2[not_qname & yes_strand & yes_tname]

    # End program if no rows
    if len(df2) == 0:
        return None, None, None

    df2["match_distance"] = np.where(
        df2["Strand_A"] == "+",  # Where: Strand is '+'
        df2["Tstart_B"] - df2["Tend_A"],  # True:  read_2_start - read_1_end
        df2["Tstart_A"] - df2["Tend_B"],  # False: read_1_start - read_2_end
    )
    # Remove reads outside of the distance parameter
    df2 = df2[(df2["match_distance"] > 0) & (df2["match_distance"] < distance)]

    # End program if no rows
    if len(df2) == 0:
        return None, None, None

    df2 = df2.drop_duplicates(
        subset=[
            "channel",
            "start_time",
            "duration",
            "next_start_time",
            "read_id",
            "next_read_id",
            "sequence_length_template",
            "next_sequence_length_template",
            "combined_length",
        ],
        keep="first",
    )
    # separate df into read groups and set index to cs to allow grouping
    cond_1 = df2["next_read_id"] == df2["read_id"].shift(-1)
    cond_2 = df2["read_id"] == df2["next_read_id"].shift(-1)
    df2["COND"] = np.where(cond_1 | cond_2, True, False)
    df2["W"] = np.where(df2["COND"].shift(1) == False, 1, 0)
    df2["cs"] = df2["W"].cumsum()

    """UNDER HERE NOT REVISED OR COMMENTED WELL"""
    # TODO: finish commenting

    if alt:
        groupby_list = ["cs", "Tname_B"]
    else:
        groupby_list = ["cs"]
    df2 = df2.set_index(groupby_list)
    df2_groupby = df2.groupby(level=groupby_list)

    # group and concatenate read ids
    df2["all_but_last"] = df2_groupby["read_id"].apply("|".join)
    df2["last_read_id"] = df2_groupby["next_read_id"].last()

    # TODO: this is the failing point, can it be cut off sooner?

    df2["cat_read_id"] = df2["all_but_last"] + "|" + df2["last_read_id"]

    # group and combine length
    df2["combined_length"] = df2_groupby["sequence_length_template"].sum()
    df2["last_length"] = df2_groupby["next_sequence_length_template"].last()
    df2["combined_length"] = df2["combined_length"] + df2["last_length"]

    # take max/min for end/start match from grouped value list
    df2["start_match"] = (
        df2_groupby[["Tstart_A", "Tstart_B", "Tend_A", "Tend_B"]]
        .transform("min")
        .min(axis=1)
    )
    df2["end_match"] = (
        df2_groupby[["Tstart_A", "Tstart_B", "Tend_A", "Tend_B"]]
        .transform("max")
        .max(axis=1)
    )

    # group and add start and end times
    df2["start_time"] = df2_groupby["start_time"].first()
    df2["next_end"] = df2_groupby["next_end"].last()

    # add the duration (time between start and end)
    df2["duration"] = df2["next_end"] - df2["start_time"]

    # format and add coordinates
    df2["stime_floor"] = np.floor(df2["start_time"]).astype("int64").astype("str")
    df2["etime_ceil"] = np.ceil(df2["next_end"]).astype("int64").astype("str")
    df2["channel"] = df2["channel"].astype("int64").astype("str")
    df2["combined_length"] = df2["combined_length"].astype("int64")
    df2["start_match"] = df2["start_match"].astype("int64").astype("str")
    df2["end_match"] = df2["end_match"].astype("int64").astype("str")
    df2["duration"] = df2["duration"].map("{:.5f}".format)
    df2["coords"] = df2["channel"] + ":" + df2["stime_floor"] + "-" + df2["etime_ceil"]

    # rename cols for export
    df2.rename(columns={"Tname_A": "target_name", "Strand_A": "strand"}, inplace=True)

    # fused_read_ids is a pd.Series of all fused reads
    fused_read_ids = pd.concat([df2["read_id"], df2["next_read_id"]])

    df2["count"] = df2_groupby.size() + 1

    # remove duplicate entries from df2
    df2 = df2.drop_duplicates(
        subset=[
            "coords",
            "channel",
            "start_time",
            "duration",
            "combined_length",
            "start_match",
            "end_match",
            "cat_read_id",
        ],
        keep="first",
    )
    fused_read_ids = fused_read_ids.unique()

    # un_fused_df contains reads that are correctly split
    un_fused_df = seq_sum_df[~seq_sum_df["read_id"].isin(fused_read_ids)].reset_index()
    # split_df is reads that have false starts (i.e 2->N)
    split_df = seq_sum_df[seq_sum_df["read_id"].isin(fused_read_ids)].reset_index()

    # TODO: CLEAN UP EXTRA SERIES FROM DFS
    return df2, un_fused_df, split_df


def die(message, status=1):
    """Print an error message and call sys.exit with the given status, terminating the process"""
    print(message, file=sys.stderr)
    sys.exit(status)


def print_args(args, label="Arguments"):
    """Print and format all arguments from the command line"""
    print(label + ":")
    dirs = dir(args)
    m = max([len(a) for a in dirs if a[0] != "_"])
    for attr in dirs:
        if attr[0] != "_":
            print("{a:<{m}}\t{b}".format(a=attr, m=m, b=getattr(args, attr)))
    print("========================================")


if __name__ == "__main__":
    sys.exit("ERROR: core is not directly executable")
