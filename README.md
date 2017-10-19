# Estimating graph properties with random walks

## Intro
Counting the number of nodes and edges of a graph is easy right? You just use a package like [NetworkX](https://networkx.github.io/documentation/networkx-1.10/reference/functions.html) and type in 'number_of_nodes(G)'.

Well. Yes, normally. But what about when you can't do that. What if the graph in question is the entire world wide web, or all videos on youtube. For one, the data may be too large and not available on any single centralised server, and secondly we might have only partial access to the graph - we can only view the currenty node that we're on (e.g. the hyperlinks of the current webpage).

So how can we estimate properties of these graphs?

We use random walks. Specifically we choose some point to start from, see how long it takes until we return to that node

## An analogy

Imagine you're stand

## Papers

