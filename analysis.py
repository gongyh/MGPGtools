import logging
from utils.taxonomy import Taxonomy
from static import Stat
from extern.odgi import Odgi
from search import Search
from describe import Describe
from tree import Tree
from exceptions import *


class Analysis(object):
    def __init__(self):
        self.logger = logging.getLogger("timestamp")
        self.warnings = logging.getLogger("warnings")

    def statistic(self, options):
        stat_args_judgment(options)
        rank_index = Taxonomy.rank_index[options.rank]
        stat = Stat(options.db, rank_index, options.name, options.outdir)
        stat.write()

    def viz(self, options):
        viz_args_judgment(options)
        odgi = Odgi(
            options.db, options.name, options.outdir, options.outName, options.threads
        )
        odgi.viz(options)

    def search(self, options):
        search = Search(options)
        search.ref_info()

    def describe(self, options):
        describe = Describe(options)
        describe.describe_gfa()

    def tree(self, options):
        tree = Tree(options)
        if options.gfa:
            tree.drawTreeWithGfa(options.gfa, options.sampleTxt, options.label)
        else:
            tree.drawTree()
