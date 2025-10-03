import os
import re
import csv
import sys
import networkx as nx
import graphviz
from collections import defaultdict

def extract_node_ordering(main_script_path, pipeline_name):
    """Extract the order of nodes from main script"""
    with open(main_script_path, 'r') as f:
        content = f.read()
    
    # Find all nextJob_*.txt files in order
    pattern = f'{pipeline_name}_nextJob_[^"\s]+\.txt'
    nodes = re.findall(pattern, content)
    return nodes

def extract_step_relationships(pipeline_dir, pipeline_name):
    """Extract thisStep and nxtStep relationships from pipeline_*.sh files"""
    relationships = defaultdict(set)
    
    # Get all pipeline_*.sh files in the pipeline directory
    pipeline_files = [f for f in os.listdir(pipeline_dir) if f.startswith(f'{pipeline_name}_') and f.endswith('.sh')]
    
    for pipeline_file in pipeline_files:
        file_path = os.path.join(pipeline_dir, pipeline_file)
        with open(file_path, 'r') as f:
            content = f.read()
            
            # Find thisStep variable
            this_step_match = re.search(r'thisStep="([^"]+)"', content)
            if this_step_match:
                this_step = this_step_match.group(1)
                
                # Find all nxtStep variables
                nxt_steps = re.findall(r'nxtStep[0-9]*="([^"]+)"', content)
                for nxt_step in nxt_steps:
                    relationships[this_step].add(nxt_step)
    
    return relationships

def create_dot_visualization(adjacency_matrix, all_steps_sorted, source_node, terminal_node, pipeline_name):
    # Create a new directed graph
    dot = graphviz.Digraph(comment=f'{pipeline_name} Pipeline DAG')
    dot.attr(rankdir='LR')  # Left to right layout
    
    # Set node attributes
    dot.attr('node', shape='box', style='rounded,filled', fillcolor='lightblue')
    
    # Add nodes
    for step in all_steps_sorted:
        # Shorten the label by removing prefix and suffix
        short_label = step.replace(f'{pipeline_name}_nextJob_', '').replace('.txt', '')
        
        # Set different colors for source and terminal nodes
        if step == source_node:
            dot.node(step, short_label, fillcolor='lightgreen')
        elif step == terminal_node:
            dot.node(step, short_label, fillcolor='lightpink')
        else:
            dot.node(step, short_label)
    
    # Add edges
    for i, step1 in enumerate(all_steps_sorted):
        for j, step2 in enumerate(all_steps_sorted):
            if adjacency_matrix[i][j] == 1:
                dot.edge(step1, step2)
    
    return dot

def main():
    # Check if required arguments are provided
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python3 build_dag.py <pipeline_dir> <pipeline_name> [output_dir]")
        print("Example: python3 build_dag.py /path/to/pipes/gabriel gabriel /path/to/output")
        sys.exit(1)
    
    pipeline_dir = sys.argv[1]
    pipeline_name = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) == 4 else os.getcwd()
    
    print(f"Pipeline directory: {pipeline_dir}")
    print(f"Pipeline name: {pipeline_name}")
    print(f"Output directory: {output_dir}")
    
    # Verify pipeline directory exists
    if not os.path.exists(pipeline_dir):
        print(f"Error: Pipeline directory '{pipeline_dir}' does not exist")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created/verified output directory: {output_dir}")
    except Exception as e:
        print(f"Error creating output directory: {e}")
        sys.exit(1)
    
    # Find the main script
    main_script_path = os.path.join(pipeline_dir, f'{pipeline_name}_main.sh')
    print(f"Looking for main script: {main_script_path}")
    if not os.path.exists(main_script_path):
        print(f"Error: Main script '{main_script_path}' does not exist")
        sys.exit(1)
    
    # Extract node ordering from main script
    print("Extracting node ordering from main script...")
    all_steps_sorted = extract_node_ordering(main_script_path, pipeline_name)
    
    if not all_steps_sorted:
        print(f"Warning: No nodes found in {main_script_path}")
        print(f"Please check if the script contains nodes in the format: {pipeline_name}_nextJob_*.txt")
        sys.exit(1)
    
    print(f"Found {len(all_steps_sorted)} nodes in main script")
    
    # Extract step relationships from pipeline_*.sh files
    print("Extracting step relationships from pipeline files...")
    step_relationships = extract_step_relationships(pipeline_dir, pipeline_name)
    print(f"Found {len(step_relationships)} step relationships")
    
    # Create adjacency matrix
    n = len(all_steps_sorted)
    adjacency_matrix = [[0] * n for _ in range(n)]
    
    # Fill adjacency matrix
    for step, next_steps in step_relationships.items():
        if step in all_steps_sorted:
            i = all_steps_sorted.index(step)
            for next_step in next_steps:
                if next_step in all_steps_sorted:
                    j = all_steps_sorted.index(next_step)
                    adjacency_matrix[i][j] = 1
    
    # Set source and terminal nodes based on pipeline name
    source_node = f'{pipeline_name}_nextJob_copyFastqs.txt'
    terminal_node = f'{pipeline_name}_nextJob_checkProjectComplete.txt'
    print(f"Source node: {source_node}")
    print(f"Terminal node: {terminal_node}")
    
    # Create visualization using Graphviz
    print("Creating visualization...")
    dot = create_dot_visualization(adjacency_matrix, all_steps_sorted, source_node, terminal_node, pipeline_name)
    visualization_path = os.path.join(output_dir, f'{pipeline_name}_dag_visualization.png')
    print(f"Saving visualization to: {visualization_path}")
    try:
        dot.render(visualization_path.replace('.png', ''), format='png', cleanup=True)
        print("Visualization saved successfully")
    except Exception as e:
        print(f"Error saving visualization: {e}")
        sys.exit(1)
    
    # Write adjacency matrix to CSV with shortened labels
    print("Writing adjacency matrix to CSV...")
    shortened_labels = [step.replace(f'{pipeline_name}_nextJob_', '').replace('.txt', '') for step in all_steps_sorted]
    output_csv = os.path.join(output_dir, f'{pipeline_name}_dag_adjacency_matrix.csv')
    print(f"Saving CSV to: {output_csv}")
    try:
        with open(output_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([''] + shortened_labels)  # Header row
            for i, step in enumerate(all_steps_sorted):
                writer.writerow([shortened_labels[i]] + adjacency_matrix[i])
        print("CSV saved successfully")
    except Exception as e:
        print(f"Error saving CSV: {e}")
        sys.exit(1)
    
    print(f"\nSummary:")
    print(f"Generated visualization: {visualization_path}")
    print(f"Generated adjacency matrix: {output_csv}")
    print(f"Found {len(all_steps_sorted)} nodes in the pipeline")

if __name__ == '__main__':
    main() 