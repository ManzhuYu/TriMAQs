# Report figures

Open Report_figures.ipynb in Jupyter. Set DATA_ROOT and choose the publication folder layout. All other paths are relative. The notebook saves high-resolution figures, captions and plotting tables under outputs/figures; it can display the same figures in the notebook.

Figure numbering follows the “Trends in Exposure to Air Pollutants from Marcellus and Utica Shale Development (TriMAQs)” report by Gernand et al. The concentration maps are Figures 13–16, the time series is Figure 17, and source-contribution maps are Figures 36–38. The data needed for these sections are in the yearly concentration and source-contribution NetCDFs. NO2/SO2 are already ppb, so the notebook does not convert them again.

Figures 5 and 6 read the final shared emission inventories directly, using tonnes per day and tonnes per month respectively. Other sections need their population/exposure, regional-emission, monitoring or sensitivity tables. They are disabled until those tables are available and reviewed. The notebook lists the required inputs for each section. Raw ACS files are not needed to plot prepared population/exposure results.

Publish code and the notebook on GitHub. Include the final figure images and captions with the ScholarSphere dataset. Notebook execution outputs can be cleared before upload; the saved figure files are independent of those embedded images.

Maps and emission graphs follow the original plotting styles. County overlays are enabled and require geopandas and the supplied boundary shapefile. Valid low concentrations are displayed at the minimum color without changing the data; missing values stay blank. Source-contribution color bars show fractions.
