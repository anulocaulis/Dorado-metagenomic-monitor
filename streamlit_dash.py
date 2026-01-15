if __name__ == "__main__":
    # Initialize the CoverageTracker
    # The cache_file stores processed data to persist across Streamlit runs
    tracker = CoverageTracker(cache_file="coverage_cache.json")

    # Example of how you might add some initial species assignments for testing
    # In a real scenario, these would be populated as BAM files are processed
    # tracker.add_species_assignment("Escherichia_coli")
    # tracker.add_species_assignment("Bacillus_subtilis")
    # tracker.add_species_assignment("Escherichia_coli")
    # tracker.add_species_assignment("Salmonella_enterica")

    # To run the Streamlit dashboard:
    # 1. Save the entire notebook content (definitions and this execution block) as a Python file, e.g., 'realtime_metagenomic_monitor.py'
    # 2. Run from your terminal: `streamlit run realtime_metagenomic_monitor.py`
    # Note: Running Streamlit directly within a Colab cell like this often doesn't display the full UI correctly.
    # It's intended to be run as a standalone app.

    # Display the Streamlit dashboard (this function contains the st.set_page_config, st.title, etc.)
    tracker.display_streamlit_dashboard()
