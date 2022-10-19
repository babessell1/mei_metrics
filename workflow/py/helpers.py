import yaml
import itertools
import os
from tabulate import tabulate

with open("config/config.yaml", 'r') as handle:
    try:
        config = yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        print(exc)


def string_to_list(stringed_list: str) -> list[str]:  # Ex. "ABC, DEF, GHI"
    try:
        items = stringed_list.replace(" ", "").split(",")
    except:
        raise SyntaxError(
            'Please set multiple values for an item in config.yaml as a string in the format: '
            '"ABC, DEF, GHI, ..."')

    return items


def get_samp_id(sample_info_filepath: str, germ=False) -> list[str]:
    with open(sample_info_filepath) as handle:
        ids = [line.split("\t")[0] for line in handle.readlines()]
    if germ:
        ids.append("germ")

    return ids


def get_barcode(sample_info_filepath: str, germ=False) -> list[str]:
    with open(sample_info_filepath) as handle:
        bars = [line.split("\t")[1] for line in handle.readlines()]
    if germ:
        bars.append("germ")

    return bars


def get_filter_type(sample_info_filepath: str, germ=False) -> list[str]:
    with open(sample_info_filepath) as handle:
        ids = [line.split("\t")[2] for line in handle.readlines()]
    if germ:
        ids.append("germline")

    return ids


def get_chromosomes(stringed_list_chroms: str) -> list[str]:
    chroms = string_to_list(stringed_list_chroms)

    return list(set(chroms))


def get_mei_type(stringed_list_meis: str) -> list[str]:
    meis = [mei.upper() for mei in string_to_list(stringed_list_meis)]
    if any(mei not in ["LINE", "ALU", "SVA", "HERVK"] for mei in meis):
        raise ValueError("Accepted MEI types are LINE, ALU, SVA, HERVK")

    return list(set(meis))


def get_filter_type_bad(stringed_list_filter: str, germ=False) -> list[str]:
    filters = string_to_list(stringed_list_filter)
    if any(filter_ not in ["raw", "hp1", "hp2", "hp_un", "hp_non"] for filter_ in filters):
        raise ValueError("Accepted filter types are hp1, hp2, un, np")
    if germ:
        filters.append("germline")

    return list(set(filters))


def exp_samp_ids(samp_ids: list[str]) -> list[str]:  # expand sample ids to generate zipped wildcards with all combos
    return samp_ids*\
        len(get_mei_type(config["MEI"]))*\
        len(get_chromosomes(config["CHROMOSOMES"]))


def exp_barcodes(barcodes: list[str]) -> list[str]:
    return barcodes*\
        len(get_mei_type(config["MEI"]))*\
        len(get_chromosomes(config["CHROMOSOMES"]))


def exp_filter_types(filter_types: list[str]) -> list[str]:
    return filter_types*\
       len(get_chromosomes(config["CHROMOSOMES"]))*\
       len(get_mei_type(config["MEI"]))


def exp_chromosomes(chromosomes: list[str], germ=False, phased=False) -> list[str]:
    info_filepath_str = config["RAW_INFO_FILEPATH"] if not phased else config["PHASED_INFO_FILEPATH"]
    return list(
        itertools.chain.from_iterable(
            itertools.repeat(x, len(get_samp_id(info_filepath_str, germ=germ))
                             ) for x in chromosomes))*len(get_mei_type(config["MEI"]))


def exp_meis(meis: list[str], germ=False, phased=False) -> list[str]:
    info_filepath_str = config["RAW_INFO_FILEPATH"] if not phased else config["PHASED_INFO_FILEPATH"]
    return list(  # ['ALU', 'LINE'] -> ['ALU', 'ALU', 'ALU', 'LINE', 'LINE', 'LINE']
        itertools.chain.from_iterable(
            itertools.repeat(x, len(get_samp_id(info_filepath_str, germ=germ))
                             ) for x in meis))*len(get_chromosomes(config["CHROMOSOMES"]))


def get_bam_dir(filt: str) -> str:
    return config["RAW_DIR"] if filt == "raw" else config["PHASED_DIR"]



