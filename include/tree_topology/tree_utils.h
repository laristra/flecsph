#pragma once

namespace flecsi {
namespace topology {

/**
 * @brief Traverse the tree and draw by levels in tikz
 **/
void
tikz_draw(const char * prefix, int level) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // Create the file
  char filename[64];
  sprintf(filename, "%s_%05d_%05d.tex", prefix, level, rank);
  std::ofstream output;
  output.open(filename);

  // Create the file header
  output << "\\documentclass{standalone}" << std::endl;
  output << "\\usepackage{tikz}" << std::endl;
  output << "\\begin{document}" << std::endl;
  output << "\\begin{tikzpicture}" << std::endl;

  // Output the tree
  std::stack<hcell_t *> stk;
  stk.push(root());

  std::vector<hcell_t *> queue;
  std::vector<hcell_t *> nqueue;
  queue.push_back(root());
  int clevel = -1;
  while(!queue.empty()) {
    for(hcell_t * e : queue) {
      hcell_t * cur = e;
      key_t nkey = cur->key();
      if(cur->is_node()) {
        point_t c = cofm_[cur->node_idx()].coordinates();
        element_t r = cofm_[cur->node_idx()].radius();
        element_t l = cofm_[cur->node_idx()].lap();
        if(clevel == level - 1) {
          output << "\\draw[blue] (" << c[0] << "," << c[1]
                 << ") circle (0.005cm);" << std::endl;
          output << "\\draw[blue!60!white] (" << c[0] << "," << c[1]
                 << ") circle (" << l + r << "cm);" << std::endl;
        }
        else {
          output << "\\draw[blue,opacity=0.2] (" << c[0] << "," << c[1]
                 << ") circle (0.005cm);" << std::endl;
          output << "\\draw[blue!60!white,opacity=0.2] (" << c[0] << "," << c[1]
                 << ") circle (" << l + r << "cm);" << std::endl;
        }
        for(int i = 0; i < nchildren_; ++i) {
          if(cur->get_child(i)) {
            key_t ckey = nkey;
            ckey.push(i);
            auto it = htable_.find(ckey);
            nqueue.push_back(&(htable_.find(ckey)->second));
          } // if
        }
      }
      else {
        point_t c = entities_[cur->entity_idx()].coordinates();
        element_t r = entities_[cur->entity_idx()].radius();
        if(clevel == level - 1) {
          output << "\\draw[red] (" << c[0] << "," << c[1]
                 << ") circle (0.005cm);" << std::endl;
          output << "\\draw[red] (" << c[0] << "," << c[1] << ") circle (" << r
                 << "cm);" << std::endl;
        }
        else {
          output << "\\draw[red,opacity=0.2] (" << c[0] << "," << c[1]
                 << ") circle (0.005cm);" << std::endl;
          output << "\\draw[red,opacity=0.2] (" << c[0] << "," << c[1]
                 << ") circle (" << r << "cm);" << std::endl;
        }
      }
    }
    queue = nqueue;
    nqueue.clear();
    clevel++;
  } // while

  // Finish the file
  output << "\\end{tikzpicture}" << std::endl;
  output << "\\end{document}" << std::endl;
  output.close();
}

/**
 * @brief      Export to a file the current tree in memory
 * This is useful for small number of particles to see the tree
 * representation
 */
void
graphviz_draw(int num) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  log_one(trace) << rank << " outputing tree file #" << num << std::endl;

  char fname[64];
  sprintf(fname, "output_graphviz_%02d_%02d.gv", rank, num);
  std::ofstream output;
  output.open(fname);
  output << "digraph G {" << std::endl << "forcelabels=true;" << std::endl;

  // Add the legend
  // output << "branch [label=\"branch\" xlabel=\"sub_entities,owner\"]"
  //       << std::endl;

  std::stack<hcell_t *> stk;
  // Get root
  stk.push(root());

  while(!stk.empty()) {
    hcell_t * cur = stk.top();
    stk.pop();
    if(cur->is_unset()) {
      output << std::oct << cur->key() << std::dec << " [label=\"" << std::oct
             << cur->key() << std::dec << "\", xlabel=\"\"];" << std::endl;
      output << std::oct << cur->key() << std::dec
             << " [shape=circle,color=black]" << std::endl;

      // Add the child to the stack and add for display
      for(size_t i = 0; i < nchildren_; ++i) {
        if(cur->get_child(i)) {
          key_t ckey = cur->key();
          ckey.push(i);
          auto it = htable_.find(ckey);
          if(it != htable_.end()) {
            stk.push(&it->second);
          }
          else {
            continue;
          }
          output << std::oct << cur->key() << "->" << it->second.key()
                 << std::dec << std::endl;
        } // if
      }
    }
    else if(cur->is_node()) {
      int sub_ent = 0;
      int idx = cur->node_idx();
      cofm_t * c = cur->is_shared() ? &shared_nodes_[idx] : &cofm_[idx];
      output << std::oct << cur->key() << std::dec << " [label=\"" << std::oct
             << cur->key() << std::dec << "\", xlabel=\"" << cur->nchildren()
             << "," << c->sub_entities() << "," << cur->owner() << "\"];"
             << std::endl;
      if(cur->is_shared()) {
        output << std::oct << cur->key() << std::dec
               << " [shape=circle,color=green]" << std::endl;
      }
      else {
        output << std::oct << cur->key() << std::dec
               << " [shape=circle,color=blue]" << std::endl;
      }

      // Add the child to the stack and add for display
      for(size_t i = 0; i < nchildren_; ++i) {
        if(cur->get_child(i)) {
          key_t ckey = cur->key();
          ckey.push(i);
          auto it = htable_.find(ckey);
          if(it != htable_.end()) {
            stk.push(&it->second);
          }
          else {
            continue;
          }
          output << std::oct << cur->key() << "->" << it->second.key()
                 << std::dec << std::endl;
        } // if
      }
    }
    else {
      output << std::oct << cur->key() << std::dec << " [label=\"" << std::oct
             << cur->key() << std::dec << "\", xlabel=\"" << cur->owner()
             << "\"];" << std::endl;
      if(cur->is_shared()) {
        output << std::oct << cur->key() << std::dec
               << " [shape=circle,color=grey]" << std::endl;
      }
      else {
        output << std::oct << cur->key() << std::dec
               << " [shape=circle,color=red]" << std::endl;
      }
    } // if
  } // while
  output << "}" << std::endl;
  output.close();
}

} // namespace topology
} // namespace flecsi