import csv
import networkx as nx
from networkx.algorithms import bipartite
import os
from collections import defaultdict
import warnings
import logging

logging.basicConfig(format='%(message)s')
logger = logging.getLogger()
logger.setLevel(logging.INFO)  # set to DEBUG for debug logging

warnings.simplefilter(action='ignore', category=FutureWarning)



# ========== Part 1 ===========================
logger.info("Part 1:")

logger.debug("    loading data")

with open(os.path.join("data", "diseases.tsv"), "r", encoding='UTF-8') as f:
    rows = csv.reader(f, delimiter="\t")
    next(rows) 
    next(rows)
    data = [(r[0], [x.strip() for x in r[2].split(",")]) for r in rows]
edges = [(x[0], y) for x in data for y in x[1]]  # (disease_id, gene_code)

logger.debug(f"    Disease-Gene data sample:\n    {edges[0:3]} \n")

logger.debug("    computing stats")
stats = {}
stats["number of genes in gene-disease DB"] = len(set([x[1] for x in edges]))
stats["number of diseases in gene-disease DB"] = len(set([x[0] for x in edges]))
stats["number of links in gene-disease DB"] = len(set(edges))

logger.debug("    computing network")
G = nx.from_edgelist(edges, create_using=None)
logger.debug(f"    Disease-Gene NW stats:\n    Nodes: {G.number_of_nodes()}\n    Edges:{G.number_of_edges()}\n")

# get the genes with the highest degree (the ones involved in the must number of diseases)
b = bipartite.degrees(G, nodes=set([x[0] for x in edges]))[0]
b = sorted(b, key=lambda x: x[1], reverse=True)
logger.info(f"    Genes with highest degree: {b[0:3]}")

# compute the Gene-Gene and Disease-Disease networks
M = bipartite.biadjacency_matrix(G, row_order=set([x[0] for x in edges]))
MT = M.transpose()
GenGen = MT.dot(M)
if GenGen.shape != (stats["number of genes in gene-disease DB"], stats["number of genes in gene-disease DB"]):
    raise ValueError('Dimensions of GenGen matrix misfit')
DisDis = M.dot(MT)
if DisDis.shape != (stats["number of diseases in gene-disease DB"], stats["number of diseases in gene-disease DB"]):
    raise ValueError('Dimensions of DisDis matrix misfit')

if not os.path.exists("results"):
    os.makedirs("results")

# save the edgelists such that it can be imported in external programs (e.g. gephy) for plotting
logger.debug("    saving results")
if not os.path.exists(os.path.join("results", "disdis.csv")) or not os.path.exists(os.path.join("results", "gengen.csv")):
    for adj, name in zip([GenGen, DisDis], ["gengen", "disdis"]):
        g = nx.from_numpy_matrix(adj)
        with open(os.path.join("results", name + ".csv"), "wb") as f:
            f.write("source;target;weight\n".encode('utf-8'))
            nx.write_weighted_edgelist(g, f, delimiter=";")

# ========== Part 2 ===========================
logger.info("Part 2:")

logger.debug("    loading data")
with open(os.path.join("data", "drugs.tsv"), "r", encoding='UTF-8') as f:
    rows = csv.reader(f, delimiter="\t")
    next(rows)
    drugs = [(r[0], r[7]) for r in rows]  # (gene_code, drug_name)

# compute the statistics
logger.debug("    computing stats")
g1 = set([x[1] for x in edges])
g2 = set([x[0] for x in drugs])
g3 = [x for x in g2 if x in g1]
g4 = set([x[1] for x in drugs])
stats["number of drugs in drug-gene DB"] = len(g4)
stats["number of genes in drug-gene DB"] = len(g2)
stats["number of links in drug-gene DB"] = len(drugs)
stats["number of genes in both DBs"] = len(g3)

logger.debug("    computing network")
gene_disease_links = defaultdict(list)
for e in edges:
    gene_disease_links[e[1]].append(e[0])

gd_edges = [(y, x[1]) for x in drugs for y in gene_disease_links[x[0]]]
G = nx.from_edgelist(gd_edges, create_using=None)
stats["number of drugs in drug-disease DB"] = len(set([x[1] for x in gd_edges]))
stats["number of diseases in drug-disease DB"] = len(set([x[0] for x in gd_edges]))
stats["number of links in drug-disease DB"] = len(set(gd_edges))

# compute disease-disease matrix coupled by drugs
M = bipartite.biadjacency_matrix(G, row_order=set([x[0] for x in gd_edges]))
MT = M.transpose()
dis_dis = M.dot(MT)

# write to csv file
logger.debug("    saving results")
if not os.path.exists(os.path.join("results", "disdis_drugs.csv")):
    g = nx.from_numpy_matrix(dis_dis)
    with open(os.path.join("results", "disdis_drugs.csv"), "wb") as f:
        f.write("source;target;weight\n".encode('utf-8'))
        nx.write_weighted_edgelist(g, f, delimiter=";")

for s in stats:
    logger.info(f"    {s}: {stats[s]}")

logger.info("Network edgelists are in ./results")
