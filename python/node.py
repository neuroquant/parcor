import networkx as nx
import argparse
import sys

def get_nodes(infile):
	nodes = []
	f = open(infile)
	first_line = True
	for line in f:
		words = line.split()
		if first_line:
			first_line = False
			continue
		nodes.append(words[0])
	return nodes

def parse_input(infile, directed, weighted):
	if directed:
		g = nx.DiGraph()
	else:
		g = nx.Graph()
	f = open(infile)
	nodes = get_nodes(infile)
	first_line = True
	for line in f:
		if first_line:
			first_line = False
			continue
		words = line.split()
		g.add_node(words[0])
		for i in range(1,len(words)):
			weight = float(words[i])
			if weight == 0:
				continue
			if weighted:
				g.add_edge(words[0], nodes[i-1], weight=weight)
			else:
				g.add_edge(words[0], nodes[i-1])
	return g, nodes

def node_measures(g, measures, directed, normalize):
	if not directed:
		if len(g) > 100:
			print(len(g))			
			print('Using approximate current flow betweenness as p > 100')
			measures['betweenness centrality'] = nx.approximate_current_flow_betweenness_centrality(g,weight='weight',normalized=normalize)
			measures['closeness centrality'] = nx.current_flow_closeness_centrality(g)
		else:	
			measures['betweenness centrality'] = nx.current_flow_betweenness_centrality(g,weight='weight',normalized=normalize)
			measures['closeness centrality'] = nx.current_flow_closeness_centrality(g)
			measures['closeness vitality'] = nx.closeness_vitality(g,weight='weight')
	else:	
		h,a = nx.hits_numpy(g,normalized=normalize)
		measures['hub pagerank'] = h
		measures['authority pagerank'] = a

def output_measures(nodes, measure_names, measures, outfile):
	if outfile is None:
		out = sys.stdout
	else:
		out = open(outfile,'w')
	out.write('\t')
	out.write('\t'.join(measure_names))
	out.write('\n')
	for n in nodes:
		out.write(n)
		for measure in measure_names:
			out.write('\t')
			out.write(str(measures[measure][n]))
		out.write('\n')

def main():
	if len(sys.argv) <= 1:
		infile = raw_input('Input file? ')
		outfile = raw_input('Output file (will be overwritten)? ')
		directed = (raw_input('Directed (d) or undirected?') == 'd')
		weighted = (raw_input('Weighted (w) or unweighted?') == 'w')
		normalize = (raw_input('Normalized (n) or unnormalized?') == 'w')
	else:
		parser = argparse.ArgumentParser(description='Process and analyze a graph\'s node properties.')
		parser.add_argument('input', metavar='infile', type=str,
									 action='store', help='input file name')
		parser.add_argument('output', metavar='outfile', type=str,
									 action='store', help='output file name')
		parser.add_argument('-d', '--directed',
									 action='store_true', dest='dir', help='toggle directed graph')
		parser.add_argument('-w', '--weighted',
									 action='store_true', dest='wts', help='toggle weighted graph')
		parser.add_argument('-n', '--normalized',
									 action='store_true', dest='norm', help='turn on normalization for measures')
		args = parser.parse_args()
		infile = args.input
		outfile = args.output
		directed = args.dir
		weighted = args.wts
		normalize = args.norm

	g, nodes = parse_input(infile, directed, weighted)
	measures = {}
	if len(g) > 100:
		measure_names = ['betweenness centrality', 'closeness centrality']	
	else:
		measure_names = ['betweenness centrality', 'closeness centrality', 'closeness vitality']		
	if directed:
		measure_names = ['hub pagerank', 'authority pagerank']
	node_measures(g, measures, directed, normalize)
	
	output_measures(nodes, measure_names, measures, outfile)

main()