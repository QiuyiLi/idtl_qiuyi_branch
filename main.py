from shared_utility import *
from species_tree import *
from gene_tree import *

import pprint
import os, shutil
import numpy as np
import skbio
import time
import getopt, sys


def build_tree_recurse(gene_tree, path):
    loss_nodes = []
    current_tree = gene_tree
    path_splited = path.split('_')
    _id = '_' + path.split('_')[-1] if len(path_splited) > 1 else ''
    parent_name_splited = gene_tree.name.split('_')
    _parent_id = '_' + parent_name_splited[-1] if len(parent_name_splited) > 1 else ''
    for i in range(len(parent_name_splited)):
        if ('l' in parent_name_splited[i] and len(parent_name_splited[i])==2):
            _parent_id = '_' + parent_name_splited[i-1] if i > 1 else ''
            break

    if (_id):
        subtree_path = os.path.join(path, 'gene_tree.txt')
        f = open(subtree_path)
        subtree = skbio.read(f, format='newick', into=skbio.tree.TreeNode)
        f.close()
        current_tree = subtree
        event_path = os.path.join(path, 'event.txt')
        f = open(event_path)
        line = f.readline()
        splited = line.strip().split(',')
        node_name = splited[0]
        distance = float(splited[1])
        event_name = '_' + splited[2][0]
        event_index = '_' + splited[3]
            
        new_dt_node = skbio.TreeNode()
        child = None
        for node in gene_tree.traverse():
            if node.name == (node_name + _parent_id):
                child = node
                break
            elif (node.name):
                if '_dl' in node.name:
                    if node.name.split('_dl')[0] == (node_name + _parent_id):
                        child = node
                        break
                elif '_tl' in node.name:
                    if node.name.split('_tl')[0] == (node_name + _parent_id):
                        child = node
                        break
                elif '_il' in node.name:
                    if node.name.split('_il')[0] == (node_name + _parent_id):
                        child = node
                        break
                elif '_sl' in node.name:
                    if node.name.split('_sl')[0] == (node_name + _parent_id):
                        child = node
                        break
        parent = child.parent
        new_dt_node.name = node_name + event_name + _id
        new_dt_node.length = child.length - distance
        new_dt_node.parent = parent
        new_dt_node.children.append(child)
        child.length = distance
        child.parent = new_dt_node
        for i in range(len(parent.children)):
            if (parent.children[i].name == child.name):
                del parent.children[i]
                break
        parent.children.append(new_dt_node)
        new_dt_node.children.append(subtree)
        subtree.parent = new_dt_node   
        for node in subtree.traverse():
            node.name = node.name + _id

    files = os.listdir(path)
    files_end_with_digit = []
    for f in files:
        if f.split('_')[-1].isdigit():
            files_end_with_digit.append(f)
    files = files_end_with_digit
    files = sorted(files,key = lambda x:int(x.split('_')[-1]))
    for f in files:
        file_path = os.path.join(path, f)
        if os.path.isdir(file_path):
            loss_nodes += build_tree_recurse(current_tree, file_path)

    files = os.listdir(path)

    for f in files:
        file_name = f.split('_')
        if (file_name and file_name[0] == 'ils'):
            _index = '_' + file_name[1]
            ils_path = os.path.join(path, f)
            file_ = open(ils_path)
            line = file_.readline()
            splited = line.split(',')
            node_name = splited[0]
            for node in current_tree.traverse():
                if node.name == (node_name + _id):
                    node.name = node_name + '_i' + _index + '_id' + _id
                    break
            file_.close()

    for f in files:
        file_name = f.split('_')
        if (file_name and file_name[0] == 's'):
            _index = '_' + file_name[1]
            s_path = os.path.join(path, f)
            file_ = open(s_path)
            line = file_.readline()
            splited = line.split(',')
            node_name = splited[0]
            for node in current_tree.traverse():
                if node.name == (node_name + _id):
                    node.name = node_name + '_s' + _index + '_id' + _id
                    break
            file_.close()

    for f in files:
        file_name = f.split('_')
        if (file_name and file_name[0] == 'loss'):
            loss_path = os.path.join(path, f)
            file_ = open(loss_path)
            line = file_.readline()
            splited = line.split(',')
            node_l_name = splited[0]
            node_l_distance = float(splited[1])
            node_index = int(splited[2])

            new_l_node = skbio.TreeNode()
            child = None
            print('start')
            for node in current_tree.traverse():
                if (node.name):
                    print('n=', node.name)
                    splited = node.name.split('_')
                if (_id == ''):
                    if (splited[0] == node_l_name):
                        child = node
                        break
                else:
                    if node.name == (node_l_name + _id):
                        child = node
                        break
                    elif (splited[0] == node_l_name and ('_' + splited[-1]) == _id):
                        child = node
                        break
                    elif '_dl' in node.name:
                        if node.name.split('_dl')[0] == (node_l_name + _id):
                            child = node
                            break
                        else:
                            splited = node.name.split('_dl')[0].split('_')
                    elif '_tl' in node.name:
                        if node.name.split('_tl')[0] == (node_l_name + _id):
                            child = node
                            break
                        else:
                            splited = node.name.split('_tl')[0].split('_')
                    elif '_il' in node.name:
                        if node.name.split('_il')[0] == (node_l_name + _id):
                            child = node
                            break
                        else:
                            splited = node.name.split('_il')[0].split('_')
                    elif '_sl' in node.name:
                        if node.name.split('_sl')[0] == (node_l_name + _id):
                            child = node
                            break
                        else:
                            splited = node.name.split('_sl')[0].split('_')
                    if (splited[0] == node_l_name and ('_' + splited[-1]) == _id):
                        child = node
                        break
            print(node_l_name + _id)
            print('end')
            parent = child.parent
            if '_s_' in parent.name:
                new_l_node.name = node_l_name + '_l' + '_id' + _id 
                if (parent.children[0].name == child.name):
                    parent.children[1].name += '_sl' + '_ind_' + parent.name.split('_s_')[1].split('_')[0] + '_' + str(node_index)
                else:
                    parent.children[0].name += '_sl' + '_ind_' + parent.name.split('_s_')[1].split('_')[0] + '_' + str(node_index)
            elif '_d_' in parent.name:
                new_l_node.name = node_l_name + '_l' + '_id' + _id 
                if (parent.children[0].name == child.name):
                    parent.children[1].name += '_dl' + '_ind_' + parent.name.split('_d_')[1].split('_')[0] + '_' + str(node_index)
                else:
                    parent.children[0].name += '_dl' + '_ind_' + parent.name.split('_d_')[1].split('_')[0] + '_' + str(node_index)
            elif '_t_' in parent.name:
                new_l_node.name = node_l_name + '_l' + '_id' + _id 
                if (parent.children[0].name == child.name):
                    parent.children[1].name += '_tl' + '_ind_' + parent.name.split('_t_')[1].split('_')[0] + '_' + str(node_index)
                else:
                    parent.children[0].name += '_tl' + '_ind_' + parent.name.split('_t_')[1].split('_')[0] + '_' + str(node_index)
            elif '_i_' in parent.name:
                new_l_node.name = node_l_name + '_l' + '_id' +  _id 
                if (parent.children[0].name == child.name):
                    parent.children[1].name += '_il' + '_ind_' + parent.name.split('_i_')[1].split('_')[0] + '_' + str(node_index)
                else:
                    parent.children[0].name += '_il' + '_ind_' + parent.name.split('_i_')[1].split('_')[0] + '_' + str(node_index)

            new_l_node.length = child.length - node_l_distance
            new_l_node.parent = parent
            new_l_node.children.append(child)
            child.length = node_l_distance
            child.parent = new_l_node
            for i in range(len(parent.children)):
                if (parent.children[i].name == child.name):
                    del parent.children[i]
                    break
            parent.children.append(new_l_node)
            loss_nodes.append(new_l_node)
            file_.close()

    return loss_nodes

def build_tree(gene_tree, path):
    return build_tree_recurse(gene_tree, path)

def cut_tree(final_tree, loss_nodes):
    final = final_tree.deepcopy()
    for node in loss_nodes:
        final.remove_deleted(lambda x: x.name == node.name)
    final.prune()
    return final

def main(input, coalescentArgs, duplicationArgs, transferArgs, lossArgs,
             hemiplasy, recombination):
    shutil.rmtree('./output')
    os.mkdir('./output')
    # os.mkdir('./output/newick_gene_subtrees')
    # os.mkdir('./output/subtrees')
    tree_path = './output/tree'
    os.mkdir(tree_path)
    
    Debug.log_file = open('./output/log.txt', 'w')
    Debug.log(header='Log created on ' + time.ctime() + '\n')
    Debug.summary_file = open('./output/summary.txt', 'w')
    Debug.summary(header='Summary created on ' + time.ctime() + '\n')
    Debug.recon_file = open('./output/recon_table.txt', 'w')
    Debug.recon_table(header='gene_id\tspecies_id\tevent_name\n')  
    # Debug.summary(header='Events:\n') 
    # Debug.summary(header='type\tspecies_node\tclade_set\n')   

    qstree = SpeciesTree(newick_path=input) 
    # qstree = SpeciesTree(newick_path='data/tree_sample2.txt')    # read newick species tree
    Debug.save_tree_nodes(nodes=qstree.nodes, path='output/species_nodes_table.txt')

    SpeciesTree.global_species_tree = qstree
    # I want fake id
    SpeciesTree.global_species_tree.post_order_fake_id()
    SpeciesTree.lambda_coal = np.random.gamma(shape=coalescentArgs['shape'], scale=coalescentArgs['scale'], size=len(qstree.leaves))

    Debug.log(header='\nspecies_tree ascii_art:\n',
                         bodies=[qstree.skbio_tree.ascii_art()])
    Debug.log(header='\nspecies_nodes:\n',
                         bodies=qstree.nodes)
    Debug.log(header='\ncoalescent:\n')
    coalescent_process, _ = qstree.coalescent(distance_above_root=10000)      # do coalescece based on the species tree
    Debug.log(header='\ncoalescent_process:\n',
                         bodies=[coalescent_process],
                         pformat=True)
    time_sequences = qstree.time_sequences(coalescent_process=coalescent_process)       # convert to time sequence structure
    Debug.log(header='\ntime_sequences:\n',
                         bodies=[time_sequences],
                         pformat=True)
    

    qgtree = GeneTree(time_sequences=time_sequences, species_tree=qstree, coalescent_process=coalescent_process)        # construct newick coalescent tree
    if ('const' not in duplicationArgs):
        GeneTree.lambda_dup = np.random.gamma(shape=duplicationArgs['shape'], scale=duplicationArgs['scale'], size=len(qgtree.leaves))
    else:
        GeneTree.lambda_dup = np.repeat(duplicationArgs['const'], len(qgtree.leaves))

    if ('const' not in lossArgs):
        GeneTree.lambda_loss = np.random.gamma(shape=lossArgs['shape'], scale=lossArgs['scale'], size=len(qgtree.leaves))
    else:
        GeneTree.lambda_loss = np.repeat(lossArgs['const'], len(qgtree.leaves))

    if ('const' not in transferArgs):
        GeneTree.lambda_trans = np.random.gamma(shape=transferArgs['shape'], scale=transferArgs['scale'], size=len(qgtree.leaves))
    else:
        GeneTree.lambda_trans = np.repeat(transferArgs['const'], len(qgtree.leaves))

    GeneTree.recombination = recombination
    GeneTree.hemiplasy = hemiplasy

    Debug.save_tree_nodes(nodes=qgtree.nodes, path='output/gene_nodes_table.txt')
    # Debug.save_output(contents=[qgtree.skbio_tree],
    #                              path='output/newick_gene_subtrees/gene_tree.txt')
    Debug.log(header='\ngene_tree ascii_art:\n',
                         bodies=[qgtree.skbio_tree.ascii_art()])
    Debug.log(header='\ngene_nodes:\n',
                         bodies=qgtree.nodes)
    Debug.log(header='\ngene_tree dlt_process:\n')
    events = qgtree.dlt_process(distance=0)     # locate the duplication points on the coalescent tree
    Debug.log(header='\ngene_tree events:\n',
                         bodies=[events],
                         pformat=True)
    Debug.log(header='\ngene_tree dt_subtree:\n')
    qgtree.dt_subtree(coalescent_process=coalescent_process, events=events, path=tree_path)        # generate duplication subtrees
    
    Debug.log(header='\nfull_events:\n',
                         bodies=qgtree.full_events,
                         pformat=True)
    # save final gene tree
    final_result = qgtree.skbio_tree.deepcopy()
    loss_nodes = build_tree(final_result, './output/tree')

    Debug.save_output(contents=[final_result,final_result.ascii_art()],
                                 path='./output/final_result.txt')
    final_result_cut = cut_tree(final_result, loss_nodes)
    Debug.save_output(contents=[final_result_cut,final_result_cut.ascii_art()],
                                 path='./output/final_result_cut.txt')
    final_result_cut_copy = final_result_cut.deepcopy()
    for node in final_result_cut_copy.traverse():
        if node.is_tip():
            name = node.name
            name = name.split('*')[0]
            name = int(name)
            name = chr(65 + name)
            node.name = name + '_' + node.name
    final_result_cut_copy_output = str(final_result_cut_copy).replace("'", '')

    Debug.save_output(contents=[final_result_cut_copy_output],path='./output/final_result_cut_newick.txt')

    if (not final_result_cut):
        print("EXCEPTION: ALL LOST")
        Debug.log_file.close()
        Debug.summary_file.close()
        return
            
    # species tree table
    # SpeciesTree.global_species_tree.post_order_fake_id()

    Debug.summary(header='\nSpecies_tree_table:\n')
    Debug.summary(header='node_id\tclade_set\n')

    for node in SpeciesTree.global_species_tree.nodes:
        clade = node.clade
        for i in range(len(clade)):
            clade[i] = SpeciesTree.global_species_tree.get_fake_id_from_real_id(clade[i])
        Debug.summary(header=str(node.node_id) + '\t'
                                + str(node.fake_node_id) + '\t'
                                + str(clade) + '\n')

    if (not final_result_cut):
        print("EXCEPTION: ALL LOST")
    else:
        # gene tree table
        qgtree.construct_final_gene_nodes(final_result_cut)
        qgtree.post_order_fake_id()
        
        Debug.summary(header='\nGene_tree_table:\n')
        Debug.summary(header='gene_node_id\tclade_set\n')
        for node in qgtree.nodes:      
            clade = node.name.split('*')[:-1]
            for i in range(len(clade)):
                clade[i] = SpeciesTree.global_species_tree.get_fake_id_from_real_id(clade[i])                     
            Debug.summary(header=str(node.node_id) + '\t'
                                + str(node.fake_node_id) + '\t'
                                + str(clade) + '\n' )
    
        Debug.summary(header='\nevent_table:\n')
        Debug.summary(header='\tgene_node_id\tclade_set\tevent\tspecies_node_id\n')
        for node in qgtree.nodes:
            if ('_d_' in node.name):
                node.events.append('D')
            elif ('_t_' in node.name):
                node.events.append('T')
            elif ('_i_' in node.name):
                node.events.append('I')
            elif ('_s_' in node.name):
                node.events.append('S')
            if ('_sl_' in node.name):
                 node.events.append('SL')
            if ('_dl_' in node.name):
                 node.events.append('DL')
            if ('_tl_' in node.name):
                 node.events.append('TL')
            if ('_il_' in node.name):
                 node.events.append('IL') 
            
            index = node.name.split('_id')[0].split('_')[-1]
            find_it = False
            
            if (node.events):
                for event in node.events:
                    if (len(event)==1):
                        for element in qgtree.full_events:
                            if (str(element['index'])==index):
                                species_fake_id = SpeciesTree.global_species_tree.get_fake_id_from_real_id(element['species_node_id'])
                                find_it = True
                                break
                    elif (len(event)==2):
                        index_1 = node.name.split('_ind_')[1].split('_')[0]
                        index_2 = node.name.split('_ind_')[1].split('_')[1]
                        find_it_1 = False
                        find_it_2 = False
                        for element in qgtree.full_events:
                            if (str(element['index'])==index_1):
                                species_fake_id = SpeciesTree.global_species_tree.get_fake_id_from_real_id(element['species_node_id'])
                                species_fake_id_1 = species_fake_id
                                find_it_1 = True
                            elif (str(element['index'])==index_2):
                                species_fake_id_2 = SpeciesTree.global_species_tree.get_fake_id_from_real_id(element['species_node_id'])
                                find_it_2 = True
                            find_it = find_it_1 * find_it_2
                            if (find_it):
                                # if (species_fake_id_1 != species_fake_id_2):
                                #     print('id1=', species_fake_id_1, 'id2=', species_fake_id_2)
                                break     

            if (find_it == False):
                species_fake_id = -1
            clade = node.name.split('*')[:-2]
            for i in range(len(clade)):
                clade[i] = SpeciesTree.global_species_tree.get_fake_id_from_real_id(clade[i])
            if (node.events):
                for event in node.events:
                    Debug.event_count[event] += 1
                    if (event != 'DL'):
                        Debug.summary(header=str(node.node_id) + '\t'
                                            + str(node.fake_node_id) + '\t'
                                            + str(clade) + '\t' 
                                            + str(event) + '\t'
                                            + str(species_fake_id) + '\n' )
                        Debug.recon_table(header=str(node.fake_node_id) + '\t'
                                            + str(species_fake_id) + '\t' 
                                            + str(event) + '\n')

    Debug.summary(header='\nfull_events:\n',
                         bodies=qgtree.full_events,
                         pformat=True)

    print('Number of events: ')
    print('S: ' + str(Debug.event_count['S']))
    print('D: ' + str(Debug.event_count['D']))
    # print('Loss: ' + str(Debug.event_count['l']))
    print('T: ' + str(Debug.event_count['T']))
    print('I: ' + str(Debug.event_count['I']))
    print('SL: ' + str(Debug.event_count['SL']))
    print('TL: ' + str(Debug.event_count['TL']))
    print('IL: ' + str(Debug.event_count['IL']))
    # print('DL: ' + str(Debug.event_count['DL']))

    # node = final_result_cut.tips()[0]
    # distance_to_root = final_result_cut.distance(node)
    # for node in final_result_cut.tips():
    #     if (distance_to_root != final_result_cut.distance(node)):
    #         print("invalid output!")
    #         break
    #         return

    Debug.log_file.close()
    Debug.summary_file.close()
    return

def default(str):
    return str + ' [Default: %default]'

def parseDistributionArgs(str):
    if str == None: return {}
    pieces = str.split(',')
    opts = {}
    for p in pieces:
        if '=' in p:
            key, val = p.split('=')
        opts[key] = float(val)
    return opts

def readCommand(argv):
    """
    Processes the command used to run HIDTLModel from the command line.
    """
    from optparse import OptionParser
    usageStr = """
    USAGE:      python main.py <options>
    EXAMPLES:   (1) python main.py
                    - runs a model
                (2) python main.py --input species_tree.txt
                OR  python main.py -i species_tree.txt
    """
    parser = OptionParser(usageStr, add_help_option=False)

    parser.add_option('--help', action='store_true', help='show this help message')   

    parser.add_option('-i', '--input', dest='input',
                      help='the path to an input file of a species tree', metavar='INPUT_FILE')
    parser.add_option('-c', '--coalescentArgs', dest='coalescentArgs',
                      help=default('the parameters of the gamma distribution for coalescent, e.g., "shape=val1,scale=val2"'),
                      default='shape=3,scale=0.1')
    parser.add_option('-d', '--duplicationArgs', dest='duplicationArgs',
                      help=default('the parameters of the gamma distribution for duplication event, e.g., "shape=val1,scale=val2"'),
                      default='shape=0.5,scale=0.1')
    parser.add_option('-t', '--transferArgs', dest='transferArgs',
                      help=default('the parameters of the gamma distribution for transfer event, e.g., "shape=val1,scale=val2"'),
                      default='shape=0.3,scale=0.1')
    parser.add_option('-l', '--lossArgs', dest='lossArgs',
                      help=default('the parameters of the gamma distribution for loss event, e.g., "shape=val1,scale=val2"'),
                      default='shape=0.5,scale=0.1') 
    parser.add_option('-h', '--hemiplasy', type='int', dest='hemiplasy',
                      help=default('hemiplasy option, 0 or 1'), metavar='HEMIPLASY',
                      default=1)      
    parser.add_option('-r', '--recombination', type='int', dest='recombination',
                      help=default('recombination option, 0 or 1'), metavar='RECOMBINATION',
                      default=1)   
    

    options, otherjunk = parser.parse_args(argv)
    if len(otherjunk) != 0:
        raise Exception('Command line input not understood: ' + str(otherjunk))
    args = dict()

    if options.help:
        parser.print_help()
        sys.exit()

    # input file (a species tree in newick format)
    if not options.input:
        parser.error('The input filename not given')
    args['input'] = options.input

    # distribution arguments
    args['coalescentArgs'] = parseDistributionArgs(options.coalescentArgs)
    args['duplicationArgs'] = parseDistributionArgs(options.duplicationArgs)
    args['transferArgs'] = parseDistributionArgs(options.transferArgs)
    args['lossArgs'] = parseDistributionArgs(options.lossArgs)

    # hemiplasy option
    if options.hemiplasy != 0 and options.hemiplasy != 1:
        parser.error('Invalid hemiplasy option: ' + str(options.hemiplasy))
    args['hemiplasy'] = True if options.hemiplasy == 1 else False

    # recombination option
    if options.recombination != 0 and options.recombination != 1:
        parser.error('Invalid recombination option: ' + str(options.recombination))
    args['recombination'] = True if options.recombination == 1 else False

    return args




# def parse_arg(argv):
#     try: 
#         opts, args = getopt.getopt(argv,'r:h:',
#         ['recombination=','hemiplasy=','help'])
#     except getopt.GetoptError:
#         print('Usage: {} -r <recombination> -h <hemiplasy>'.format(sys.argv[0]))
#         sys.exit()
#     if(opts):
#         for opt, arg in opts:
#             if opt in ('-r', '--recombination'):
#                 recombination = arg
#             elif opt in ('-h', '--hemiplasy'):
#                 hemiplasy = arg
#             elif opt in ('--help'):
#                 print('Usage: {} -r <recombination> -h <hemiplasy>'.format(sys.argv[0]))
#                 sys.exit()
#     else:
#         recombination = 1
#         hemiplasy = 1

#     return recombination, hemiplasy


if __name__ == '__main__':
    # recombination, hemiplasy = parse_arg(sys.argv[1:])
    # options = {
    #     'recombination': int(recombination),
    #     'hemiplasy': int(hemiplasy)
    # }

    args = readCommand(sys.argv[1:])
    main(**args)
