from Bio import Entrez, SeqIO

Entrez.email = "A.N.Other@example.com"

taxonomy = "Apis mellifera"


def get_scientific_name(taxonomy):
    handle = Entrez.efetch(db="Taxonomy", id=taxonomy)
    record = Entrez.read(handle)
    return record[0]["ScientificName"]

handle = Entrez.esearch(db="taxonomy", term=f"{taxonomy}[next level]", retmax=9999)
record = Entrez.read(handle)
children = []
print(f"{len(record['IdList'])} children found")
for i in record["IdList"]:
    children.append(get_scientific_name(i))
print(children)


