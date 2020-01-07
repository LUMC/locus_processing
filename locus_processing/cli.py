from os import listdir
from os.path import join, isfile, basename
import warnings

import click
import yaml

from . import load_locus_yaml, dump_locus


def _validate_inputs(ctx, param, value):
    """Mutually exclusive group for -I and -D
    Returns generator of yaml paths for -D
    """
    if param.name == "input_directory":
        idir = value
        ifile = ctx.params.get("input")
    elif param.name == "input":
        ifile = value
        idir = ctx.params.get("input_directory")
    else:
        ifile = None
        idir = None

    if (ifile is None and idir is None) or \
            (ifile is not None and idir is not None):
        raise click.BadArgumentUsage("Either --input OR "
                                     "--input-directory must be set")

    if idir is not None:
        candidates = [join(idir, x) for x in listdir(idir)]
        return filter(lambda x: isfile(x) and
                      (x.endswith('yaml') or x.endswith('yml')),
                      candidates)

    return value


def _validate_output(ctx, param, value):
    if param.name == "output_directory":
        idir = ctx.params.get("input_directory")
        if idir is None and value is not None:
            raise click.BadArgumentUsage("-O must be combined with -D")
        elif idir is not None and value is None:
            raise click.BadArgumentUsage("-O must be combined with -D")
    return value


def _single_locus_to_bed(path, prefix):
    bed_fmt = "{chrom}\t{start}\t{end}\t{name}"
    loc = load_locus_yaml(path)
    if prefix is not None:
        name = "{0}_{1}".format(prefix, loc.name)
    else:
        name = loc.name
    bed = bed_fmt.format(
        chrom=loc.chromosome.name,
        start=loc.coordinates.start,
        end=loc.coordinates.end,
        name=name
    )
    return bed


@click.command(short_help="Write regions in locus files to bed")
@click.option("--input", "-I", type=click.Path(exists=True),
              help="Path to input locus file")
@click.option("--input-directory", "-D",
              type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help="Path to directory containing locus files",
              callback=_validate_inputs)
@click.option("--prefix", "-p",
              type=click.STRING, help="Prefix to region names")
def locus_to_bed(input=None, input_directory=None, prefix=None):
    if input:
        print(_single_locus_to_bed(input, prefix))
    else:
        for y in input_directory:
            print(_single_locus_to_bed(y, prefix))


def _validate_single_locus(path):
    click.echo("{0}: ".format(path), nl=False)
    try:
        _ = load_locus_yaml(path)
    except ValueError as error:
        return_str = "ERROR, {r}".format(r=str(error))
        click.secho(return_str, fg="red")
        return 1
    else:
        return_str = "PASS"
        click.secho(return_str, fg="green")
        return 0


@click.command(short_help="Validate locus definiton files")
@click.option("--input", "-I", type=click.Path(exists=True),
              help="Path to input locus file")
@click.option("--input-directory", "-D",
              type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help="Path to directory containing locus files",
              callback=_validate_inputs)
def validate_locus(input=None, input_directory=None):
    if input is not None:
        exit(_validate_single_locus(input))
    else:
        ret = 0
        for y in input_directory:
            ret += _validate_single_locus(y)
        exit(ret)


def _complete_single_locus(fpath):
    l = load_locus_yaml(fpath)
    l.apply_hgvs_descriptions(warn_on_error=True)
    return dump_locus(l)


@click.command(short_help="Attempt to complete HGVS descriptions "
                          "within locus definitions")
@click.option("--input", "-I", type=click.Path(exists=True),
              help="Path to input locus file")
@click.option("--input-directory", "-D",
              type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help="Path to directory containing locus files",
              callback=_validate_inputs)
@click.option("--output-directory", "-O",
              type=click.Path(exists=True, file_okay=False, dir_okay=True,
                              writable=True),
              help="Path to output directory", callback=_validate_output)
@click.option("--suppress-warnings", type=click.BOOL,
              help="Suppress warnings")
def complete_locus(input=None, input_directory=None, output_directory=None,
                   suppress_warnings=False):
    if suppress_warnings:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _complete_locus(input, input_directory, output_directory)
    else:
        _complete_locus(input, input_directory, output_directory)


def _complete_locus(input=None, input_directory=None, output_directory=None):
    if input is not None:
        l = _complete_single_locus(input)
        print(yaml.dump(l, default_flow_style=False))
    else:
        for p in input_directory:
            opath = join(output_directory, basename(p))
            l = _complete_single_locus(p)
            with open(opath, "w") as ohandle:
                yaml.dump(l, ohandle, default_flow_style=False)
