import warnings
from typing import List

from .lookups import fetch_rsid, fetch_sequence, position_converter,\
    name_checker
from .utils import Result


class Chromosome(object):
    def __init__(self, name: str, accession: str):
        self.name = name
        self.accession = accession

    def __repr__(self):
        return "<Chromosome(name={n}, accession={a})>".format(
            n=self.name,
            a=self.accession
        )


class Coordinates(object):
    def __init__(self, start: int, end: int):
        self.start = start
        self.end = end

    def __repr__(self):
        return "<Coordinates(start={s}, end={e})>".format(
            s=self.start,
            e=self.end
        )

    def __len__(self):
        return self.end - self.start


class Snp(object):
    def __init__(self, id: str, g_notation: str, alt_notation: str,
                 ref_g_notation: str, c_notation: str, p_notation: str,
                 description: str, tags: List[str]):
        self.id = id
        self.g_notation = g_notation
        self.alt_notation = alt_notation
        self.ref_g_notation = ref_g_notation
        self.c_notation = c_notation
        self.p_notation = p_notation
        self.description = description
        self.tags = tags

        self.__rs_id_lookup = None

    def __repr__(self):
        return "<Snp(id={i}, g_notation={g})>".format(
            i=self.id,
            g=self.g_notation
        )

    @property
    def __get_rs_id_lookup(self) -> dict:
        if self.__rs_id_lookup is None:
            self.__rs_id_lookup = fetch_rsid(self.alt_notation)
        return self.__rs_id_lookup

    @property
    def minor_allele(self) -> str:
        return self.__get_rs_id_lookup.get("minor_allele", "")

    @property
    def major_allele(self) -> str:
        return self.__get_rs_id_lookup.get("ancestral_allele", "")

    @property
    def maf(self) -> float:
        return self.__get_rs_id_lookup.get("MAF", 0.0)

    @property
    def synonyms(self) -> List[str]:
        return self.__get_rs_id_lookup.get("synonyms", [])

    def apply_hgvs(self, build: str, chromosome: str,
                   overwrite: bool, transcript: str = None) -> Result:
        """
        Apply hgvs to the c_notation and p_notation fields.

        This will require several lookups to mutalyzer.
        Returns a Result object with any values, errors and warnings.

        If result.is_successful, results will have been applied to the object.
        the result will contain a warning of c_notation and p_notation have
        already been set, and overwrite=True
        :param build: genome build
        :param chromosome: chromosome accession nr
        :param overwrite: overwrite old c_notation and p_notation
        :param transcript: preferred transcript to select
        :return:
        """
        res = Result(False, values={}, errors=[], warnings=[])
        if not overwrite and self.c_notation is not None:
            res.errors.append("overwrite is disabled and c_notation is already set")
        if not overwrite and self.p_notation is not None:
            res.errors.append("overwrite is disabled and p_notation is already set")
        if len(res.errors) > 0:
            return res

        if self.c_notation is not None:
            res.warnings.append("c_notation was already set")
        if self.p_notation is not None:
            res.warnings.append("p_notation was already set")

        variant = "{c}:g.{v}".format(c=chromosome, v=self.g_notation)
        position_response = position_converter(build, variant)
        if not position_response.ok:
            res.errors.append("Could not get transcript "
                              "position because {0}".format(
                position_response.json().get("faultstring"))
            )
            return res
        transcripts = position_response.json()
        if transcripts is None or len(transcripts) == 0:
            res.errors.append("No transcripts found")
            return res

        transcript_names = list(map(lambda x: x.split(":")[0], transcripts))
        if transcript is not None and transcript in transcript_names:
            tr = transcripts[transcript_names.index(transcript)]
        elif transcript is not None:
            res.warnings.append("Could not find preferred transcript")
            tr = transcripts[0]
        else:
            if len(transcripts) > 1:
                res.warnings.append("More than one transcript found. "
                                    "Will select first transcript")

            tr = transcripts[0]
        res.values['c_notation'] = tr

        protein_response = name_checker(tr)
        protein_json = protein_response.json()
        if not protein_response.ok:
            error = protein_json.get("messages")[0].get("message")
            res.errors.append("Could not get protein "
                              "description because {0}".format(error))
            return res

        proteins = protein_json.get("proteinDescriptions")
        if len(proteins) > 1:
            res.warnings.append("More than one protein description")
            protein = proteins[0]
        elif len(proteins) == 0:
            res.warnings.append("No protein descriptions")
            protein = None
        else:
            protein = proteins[0]

        res.values['p_notation'] = protein

        self.c_notation = tr
        self.p_notation = protein
        res.success = True
        return res


class Haplotype(object):
    def __init__(self, name: str, type: str, snps: List[str]):
        self.name = name
        self.type = type
        self.snps = snps

    def __repr__(self):
        fmt = "<Haplotype(name={n}, type={t})>".format(
            n=self.name,
            t=self.type
        )
        return fmt


class Locus(object):
    def __init__(self, version: str, name: str, reference: str,
                 chromosome: Chromosome, coordinates: Coordinates,
                 transcript: str, snps: List[Snp],
                 haplotypes: List[Haplotype]):
        self.version = version
        self.name = name
        self.reference = reference
        self.chromosome = chromosome
        self.coordinates = coordinates
        self.transcript = transcript
        self.snps = snps
        self.haplotypes = haplotypes

        self.__sequence = None

    def __repr__(self):
        return "<Locus(name={n}, version={v})>".format(
            n=self.name,
            v=self.version
        )

    @property
    def sequence(self):
        if self.__sequence is None:
            self.__sequence = fetch_sequence(self.reference,
                                             self.chromosome.name,
                                             self.coordinates.start,
                                             self.coordinates.end)
        return self.__sequence

    def get_snp(self, id: str) -> Snp:
        """Get snp by id"""
        return next((s for s in self.snps if s.id == id), None)

    def get_haplotype_snps(self, haplotype: Haplotype) -> List[Snp]:
        """Get snps for specific haplotype"""
        return [s for s in self.snps if s.id in haplotype.snps]

    def apply_hgvs_descriptions(self, overwrite: bool = False,
                                warn_on_error: bool = False):
        """Apply hgvs descriptions for snps"""
        for snp in self.snps:
            res = snp.apply_hgvs(self.reference, self.chromosome.accession,
                                 overwrite, self.transcript)
            for w in res.warnings:
                warnings.warn("{0}: {1}".format(snp, w))
            err_concat = ",".join(res.errors)
            if len(res.errors) > 0:
                if warn_on_error:
                    warnings.warn("ERROR at {0}: {1}".format(snp, err_concat))
                else:
                    raise ValueError("ERROR at {0}: {1}".format(
                        snp, err_concat)
                    )
