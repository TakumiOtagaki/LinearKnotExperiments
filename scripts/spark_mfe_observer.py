#!/usr/bin/env python3
"""Generate a diagnostic translation unit in local scratch, never in the clone.

The official Spark binary remains the MFE reference. The observer logs traceback
node energies and can exclude added pairs outside a target extension. It does
not force target pairs: a masked result is an evaluation only when the returned
pair set equals the complete target. Unmasked output must match the official CLI.
"""
from pathlib import Path
import argparse
import re

HEADER = r'''
#include <cstdlib>
struct ObserverNode { std::string state; int i,j,energy,parent,children=0; };
static std::vector<ObserverNode> observer_nodes;
static std::vector<int> observer_stack;
static std::vector<int> observer_target;
static bool observer_allowed(int i,int j) {
    return observer_target.empty() || observer_target.at(i)==j;
}
struct ObserverScope {
    int id;
    ObserverScope(const char *state,int i,int j,int energy) {
        id=observer_nodes.size();
        int parent=observer_stack.empty()?-1:observer_stack.back();
        observer_nodes.push_back({state,i,j,energy,parent,0});
        if(parent>=0) observer_nodes[parent].children+=energy;
        observer_stack.push_back(id);
    }
    ~ObserverScope() { observer_stack.pop_back(); }
};
static void observer_load_target(int n) {
    const char *target=std::getenv("SPARK_OBSERVER_TARGET");
    if(!target) return;
    std::string s(target);
    if(s.size()!=static_cast<size_t>(n)) throw std::runtime_error("target length");
    observer_target.resize(n+1,0);
    std::vector<int> pending;
    for(int i=1;i<=n;++i) {
        if(s[i-1]=='(' || s[i-1]=='[') pending.push_back(i);
        else if(s[i-1]==')' || s[i-1]==']') {
            if(pending.empty()) throw std::runtime_error("target brackets");
            int j=pending.back();pending.pop_back();
            observer_target[i]=j;observer_target[j]=i;
        } else if(s[i-1]!='.') throw std::runtime_error("target character");
    }
    if(!pending.empty()) throw std::runtime_error("target brackets");
}
static void observer_dump(const vrna_param_t *p) {
    std::cerr<<"metadata\t"<<p->temperature<<'\t'<<p->model_details.dangles
             <<'\t'<<p->MLclosing<<'\t'<<p->MLintern[1]<<'\t'<<p->MLbase<<'\n';
    for(size_t id=0;id<observer_nodes.size();++id) {
        const auto &v=observer_nodes[id];
        std::cerr<<"loop\t"<<id<<'\t'<<v.parent<<'\t'<<v.state<<'\t'<<v.i<<'\t'<<v.j
                 <<'\t'<<v.energy-v.children<<'\t'<<v.energy<<'\n';
    }
}
'''


def generate(source: Path, destination: Path, trace_repair=False):
    text = (source / 'src/Spark.cc').read_text()
    text = text.replace('static bool pseudoknot = false;', HEADER + '\nstatic bool pseudoknot = false;', 1)
    pattern = r'(void trace_(\w+)\([^;{}]+\) \{)'
    states = []

    def observe(match):
        state = match[2]
        energy = 'spark.W_[j]' if state == 'W' else 'spark.WM2_[j]' if state == 'WM2' else 'e'
        right = 'ip' if state == 'BE' else 'j'
        states.append(state)
        return match[1] + f'\n    ObserverScope observed("{state}",i,{right},{energy});'

    text = re.sub(pattern, observe, text)
    assert len(states) == 13, states
    v_guard = 'if (ptype_closing > 0 && !restricted && evaluate && pkonly) {'
    assert text.count(v_guard) == 1
    text = text.replace(v_guard, 'if (ptype_closing > 0 && !restricted && evaluate && pkonly && (paired || observer_allowed(i,j))) {')
    vp_guard = 'if (weakly_closed_ij || tree.tree[i].pair >= -1 || tree.tree[j].pair >= -1 || ptype_closing == 0) {'
    assert text.count(vp_guard) == 1
    text = text.replace(vp_guard, 'if (weakly_closed_ij || tree.tree[i].pair >= -1 || tree.tree[j].pair >= -1 || ptype_closing == 0 || !observer_allowed(i,j)) {')
    text = text.replace('Spark spark(seq, !args_info.noGC_given, restricted);',
                        'observer_load_target(n);\n    Spark spark(seq, !args_info.noGC_given, restricted);', 1)
    text = text.replace('std::string structure = trace_back(spark, tree, mark_candidates);',
                        'std::string structure = trace_back(spark, tree, mark_candidates);\n    observer_dump(spark.params_);', 1)
    # Log forward DP alternatives, which remain usable when traceback fails.
    begin = text.index('energy_t compute_BE(')
    end = text.index('energy_t compute_WMBP(', begin)
    part = text[begin:end]
    assert part.count('return val;') == 1
    part = part.replace('return val;', r'''
    std::cerr<<"BE_candidates\t"<<i<<'\t'<<ip<<'\t'<<lp<<'\t'<<l
             <<'\t'<<ILoopE(spark.S_,spark.S1_,spark.params_,ptype_closing_iip,i,ip,lp,l)
             <<'\t'<<BE_energy<<'\t'<<m1<<'\t'<<m2<<'\t'<<m3<<'\t'<<m4<<'\t'<<m5<<'\t'<<val<<'\n';
    return val;''')
    text = text[:begin] + part + text[end:]
    if trace_repair:
        # Hypothesis check only: forward DP and energy computation are unchanged.
        # Keep this binary separate from both the official CLI and the observer.
        text = text.replace('i >= spark.n_ || j >= spark.n_', 'i > spark.n_ || j > spark.n_', 1)
        begin = re.search(r'void trace_BE\([^;{}]+\) \{', text).start()
        end = text.index('const std::string &trace_back', begin)
        part = text[begin:end]
        assert part.count('spark.WIP_Bp[j - 1]') == 1
        assert part.count('spark.WIP_Bbp[j - 1]') == 1
        part = part.replace('spark.WIP_Bp[j - 1]', 'spark.WIP_[j - 1]')
        part = part.replace('spark.WIP_Bbp[j - 1]', 'spark.WIP_[j - 1]')
        text = text[:begin] + part + text[end:]
    destination.write_text(text)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--trace-repair', action='store_true', help='diagnostic hypothesis check; not an upstream fix')
    args = parser.parse_args()
    if str(args.destination.resolve()).startswith('/gs/bs/'):
        parser.error('generated translation unit must be in node-local scratch')
    generate(args.source, args.destination, args.trace_repair)
