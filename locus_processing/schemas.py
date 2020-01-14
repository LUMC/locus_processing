from marshmallow import Schema, fields, post_load, validates_schema, post_dump
from marshmallow.exceptions import ValidationError

from .models import Chromosome, Coordinates, Snp, Haplotype, Locus


class ChromosomeSchema(Schema):
    name = fields.Str()
    accession = fields.Str()

    @post_load
    def make_chromosome(self, data):
        return Chromosome(**data)


class CoordinatesSchema(Schema):
    start = fields.Int()
    end = fields.Int()

    @post_load
    def make_coordinate(self, data):
        return Coordinates(**data)

    @validates_schema
    def validate_location(self, data):
        if data['end'] < data['start']:
            raise ValidationError("End cannot be before start")


class SnpSchema(Schema):
    id = fields.Str()
    g_notation = fields.Str()
    alt_notation = fields.Str()
    ref_g_notation = fields.Str(required=False, missing=None)
    c_notation = fields.Str(allow_none=True)
    p_notation = fields.Str(allow_none=True)
    description = fields.Str(allow_none=True)
    tags = fields.List(fields.Str(), allow_none=True)

    @post_load
    def make_snp(self, data):
        return Snp(**data)

    @post_dump
    def remove_ref_g_if_none(self, data):
        if data.get("ref_g_notation") is None:
            return {k: v for k, v in data.items() if k != "ref_g_notation"}
        return data


class HaplotypeSchema(Schema):
    name = fields.Str(allow_none=True)
    type = fields.Str(allow_none=True)
    snps = fields.List(fields.Str(), allow_none=True)

    @post_load
    def make_haplotype(self, data):
        return Haplotype(**data)


class LocusSchema(Schema):
    version = fields.Str()
    name = fields.Str()
    reference = fields.Str()
    chromosome = fields.Nested(ChromosomeSchema)
    coordinates = fields.Nested(CoordinatesSchema)
    transcript = fields.Str()
    snps = fields.Nested(SnpSchema, many=True)
    haplotypes = fields.Nested(HaplotypeSchema, many=True)

    @post_load
    def make_locus(self, data):
        return Locus(**data)

    @validates_schema
    def validate_haplotypes(self, data):
        snp_ids = [x.id for x in data.get('snps', [])]

        # Make sure there are only unique snp id's
        if len(snp_ids) != len(set(snp_ids)):
            raise ValidationError("snp id's are not unique")

        hap_types = []
        hap_snps = []
        first = True
        for hap in data['haplotypes']:
            hap_types.append(hap.type)
            hap_snps.extend(hap.snps)
            for snp in hap.snps:
                if snp not in snp_ids:
                    raise ValidationError("Haplotype {0} has unknown SNP".format(hap))

            if not first and hap.snps == []:
                raise ValidationError("Haplotype without variants is not first entry")

            first = False

        # Make sure there are only unique hap id's
        if len(hap_types) != len(set(hap_types)):
            raise ValidationError("hap types's are not unique")

        # Make sure that all variants end up in a haplotype
        if set(hap_snps) != set(snp_ids):
            raise ValidationError("Some SNPs are not used in haplotypes")
