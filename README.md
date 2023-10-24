# Space-Time-Kriging

This report summarizes the outcomes of the GIMA internship at the Centre for Geographical Analysis (CGA), a research institution associated to the University of Stellenbosch in South Africa. The CGA specializes in GIS and strives to advance spatial information management for the benefit of the Southern African community.

Climate change is becoming a growing concern in South Africa, specifically in the southwestern Cape. Over the last five decades, the mean annual temperature increased significantly while also extreme rainfall events happened more frequently. As this trend is expected to continue, adaptation strategies will be essential to mitigate the risks and safeguard agricultural production in South Africa's future. To achieve this, the provision of accurate data of local climate surfaces is crucial.

This internship project aimed to research and develop spatially interpolated climate surfaces from weather station data in support of the TerraClim agricultural decision support system (ADSS) in the southwestern Cape of South Africa. The internship was divided into two projects, each having at a different purpose within researching and developing climate surfaces. 

![FlowChart_Methodology_Transparant drawio](https://github.com/renswvw/Space-Time-Kriging/assets/94464752/52306308-621f-40f3-bf84-a5526dbbe10b)
![TerraClimRegion_WC_R](https://github.com/renswvw/Space-Time-Kriging/assets/94464752/2010035e-0f70-4e49-a810-7d087cfa7f6b)
The first project provided a theoretical foundation for exploring various climate interpolation techniques and their strengths and limitations. The knowledge obtained from this project played a vital role in selecting the appropriate interpolation technique for the second project. The foremost result was the identification of space-time kriging (STK) as the most suitable technique for the TerraClim ADSS. STK has ability to incorporate both spatial and temporal dimensions into the interpolation process, and this enables to provision of accurate real-time data for precision agriculture. Implementing this technique could facilitate more informed decision-making and improve climate monitoring in the complex terrains of the southwestern Cape.

The second project focused on the research and development of an advanced interpolation technique. This involved a literature review into STK and the development of a R script capable of implementing STK within the TerraClim ADSS. The script its adaptability and scalability for various datasets and parameters makes it useful for multiple applications. While this script was a significant achievement, it also encountered several practical challenges, including limited computational resources and the absence of a validation process.

![WeatherStations_WC_R](https://github.com/renswvw/Space-Time-Kriging/assets/94464752/38ff84ce-e185-4e18-913a-9a55a91ea7c3)
![Stations_R](https://github.com/renswvw/Space-Time-Kriging/assets/94464752/f7ae48f2-79fb-437e-a43f-535aca4d5056)
![krigeST](https://github.com/renswvw/Space-Time-Kriging/assets/94464752/6d50353e-d7e9-4d77-805d-3438bdad6815)

In short, this project provided both a theoretical and practical foundation. Despite encountered challenges, the project highlighted the potential of STK to enhance the accuracy of spatially interpolated climate surfaces within the TerraClim ADSS. This provided a solid foundation for future projects to enhance the accuracy of TerraClim ADSS, thereby improving the capacity to support local farmers in making informed decisions regarding precision agriculture and natural resource management in the southwestern Cape of South Africa.	 
