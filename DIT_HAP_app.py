import streamlit as st

plot_page = st.Page("pages/depletion_data.py", title="Curve plot", icon=":material/timeline:")
feature_space_page = st.Page("pages/feature_space.py", title="Feature space", icon=":material/scatter_plot:")
enrichment_page = st.Page("pages/enrichment_analysis.py", title="Enrichment analysis", icon=":material/search_insights:")

pg = st.navigation(
        {
            "Visualization": [plot_page, feature_space_page],
            "Analysis": [enrichment_page],
        }
    )
pg.run()


