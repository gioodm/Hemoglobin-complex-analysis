#!/usr/bin/python3

import sys
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

def parse_mitab(mitabfile, organism):
	list_int=[]
	f=open(mitabfile)
	for line in f:
		v=line.split('\t')
		if (v[0].find('uniprotkb:') == 0 and v[1].find('uniprotkb:') == 0):
			if (v[9].find(organism) > -1 and v[10].find(organism) > -1):
				int1 = v[0].replace('uniprotkb:','')
				int2 = v[1].replace('uniprotkb:','')
				texp = v[6].replace('psi-mi:','')
				cexp = v[11].replace('psi-mi:','')
				tax1 = v[9].replace('taxid:','')
				tax2 = v[10].replace('taxid:','')
				pint = (int1,int2,texp,tax1,tax2,cexp)
				list_int.append(pint)
	return list_int

def unique_interaction(filter_list):
	inter = {}
	for pint in filter_list:
		k = [pint[0],pint[1]]
		k.sort()
		inter[tuple(k)] = inter.get(tuple(k),0)+1
	return inter

def parse_regulondb(dic, l):
	g = nx.Graph()

	for k1, k2 in list(dic.keys()):
		if k1 in l:
			g.add_node(k1, color = '#E83813', node_size = 400, with_labels = True)
		else: 
			g.add_node(k1, color = '#45B39D', node_size = 15, with_labels = False)
		if k2 in l: 
			g.add_node(k2, color = '#E83813', node_size = 400, with_labels = True)
		else:
			g.add_node(k2, color = '#45B39D', node_size = 15, with_labels = False) ##labels not working properly
		g.add_edge(k1, k2, color = '#CCD1D1')

	return g

def show_graph(g,l):
	pos = nx.spring_layout(g)
	lc, ls, ll = [], [], {}
	for i in g.nodes():
		lc.append(g.nodes()[i]['color'])
		ls.append(g.nodes()[i]['node_size'])
		if i in l:
			ll[i] = i 
	ec = [g[u][v]['color'] for u,v in g.edges()]
	nx.draw(g, node_color = lc, node_size = ls, with_labels = False, edge_color = ec, pos = pos) #inform that nodes have different colours
	nx.draw_networkx_labels(g, pos, ll, font_weight = 'bold', font_family = 'verdana', font_size = 10)
	plt.show()

		
			
if __name__ == '__main__':
	mitabfile = sys.argv[1]
	organism = sys.argv[2]
	subnet = parse_mitab(mitabfile, organism)
	unique = unique_interaction(subnet)

	HbA = 'P69905'
	HbB = 'P68871'

	G = parse_regulondb(unique, [HbA, HbB])
	gA2 = nx.single_source_shortest_path_length(G, source = HbA, cutoff = 2)
	gB2 = nx.single_source_shortest_path_length(G, source = HbB, cutoff = 2)
	gAA2 = set(gA2.keys())
	gBB2 = set(gB2.keys())
	new_nodes2 = gAA2.union(gBB2)
	new_nodes2 = list(new_nodes2)

	g2 = G.subgraph(new_nodes2)
	print("Number of network nodes: " + str(len(g2.nodes())))
	print("Number of network edges: " + str(len(g2.edges())))

	dg = g2.degree()
	print("Degree alpha subunit: " + str(dg[HbA]))
	print("Degree beta subunit: " + str(dg[HbB]))
	
	cluster = nx.clustering(g2)
	print("Transitivity alpha subunit: " + str(cluster[HbA]))
	print("Transitivity beta subunit: " + str(cluster[HbB]))

	bc = nx.betweenness_centrality(g2)
	print("Betweenness centrality alpha subunit: " + str(bc[HbA]))
	print("Betweenness centrality beta subunit: " + str(bc[HbB]))

	show = show_graph(g2, [HbA, HbB])

	gA1 = nx.single_source_shortest_path_length(G, source = HbA, cutoff = 1)
	gB1 = nx.single_source_shortest_path_length(G, source = HbB, cutoff = 1)
	gAA1 = set(gA1.keys())
	gBB1 = set(gB1.keys())
	new_nodes1 = gAA1.union(gBB1)
	new_nodes1 = list(new_nodes1)

	print("Degree:")
	list_dg = []
	for node in new_nodes1:
		list_dg.append((node, dg[node]))

	list_dg.sort(key = lambda x: x[1], reverse = True)
	for n in range(len(list_dg)):
		print(str(n+1) + " " + str(list_dg[n][0]) + " " + str(list_dg[n][1]))

	print("Transitivity:")
	list_tr = []
	for node in new_nodes1:
		list_tr.append((node, cluster[node]))

	list_tr.sort(key = lambda x: x[1], reverse = True)
	for n in range(len(list_tr)):
		print(str(n+1) + " " + str(list_tr[n][0]) + " " + str(round(list_tr[n][1], 4)))

	print("Betweenness centrality:")
	list_bc = []
	for node in new_nodes1:
		list_bc.append((node, bc[node]))

	list_bc.sort(key = lambda x: x[1], reverse = True)
	for n in range(len(list_bc)):
		print(str(n+1) + " " + str(list_bc[n][0]) + " " + str(round(list_bc[n][1], 4)))




