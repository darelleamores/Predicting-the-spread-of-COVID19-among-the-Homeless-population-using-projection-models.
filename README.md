# Predicting the spread of COVID19 among the Homeless population using projection models: A focus on the unsheltered of San Bernardino County 

# Abstract

Background: The San Bernardino homeless population being affected by the COVID-19 epidemic. 

Purpose: To using a simulation using County Data to help predict the number of cases to help determine and prepare on how to allocate Public health resources in light of the COVID-19 pandemic.

Method: Rstudio using EpiModel, Public health County data and creating simulations based on COVID-19 Public health interventions.

Result: Showed that there are increased infections and cases if there are no PH interventions involved in the homeless community in SB County. Public interventions decrease adverse outcomes such as infected individuals and case fatalities.

# Public Health Issue 

The homeless population is exposed to many illness and diseases that might spread in homeless encampments. With the lack of access and the COVID-19 pandemic, it’s important to run a simulation in order to help report, allocate public health resources in order to decrease issues such as increase infection rates and hospital overflow. 

# Background of Organizations 

San Bernardino County Public Health Department 
The Public Health Department of San Bernardino (SBD) is the public health department of San Bernardino County. San Bernardino County consists of 47 cities including San Bernardino, Ontario, Redlands, Loma Linda, and Chino Hills. The Communicable Disease Section (CDS) of SBDPH “…has the responsibility of monitoring and controlling communicable disease in the County (6).” They work local organizations in the county of San Bernardino such as hospitals. CDS works to contain and track diseases by educating, tracing, and monitoring disease trends.

# Literature Review 

Covid-19
As of June 2020, there has been an estimate over a million cases and about 100,000 deaths in the United States from the Corona virus or otherwise known as COVID-19 (CDC, 2020). The origin of COVID-19 is not for certain, but studies have may have indicated that the virus mutated between bats and humans and the first cases of COVID-19 was linked to a marketplace in Wuhan, China (Shereen, Khan, Kazmi, Bashir, & Siddique, 2020). The virus is a positive strain RNA virus that is encapsulated with spikes that resembles a crown, hence the name, Corona virus(Mohammadi, Meskini, & do Nascimento Pinto, 2020). The crown on the virus serves as a receptor for attachment in the host cells that it infects (Astuti & Ysrafil, 2020). The virus then causes upper respiratory tract infections similar to a common cold and progression of the virus in the body can cause pulmonary tissue damage in the lungs, decreased immune function due to exhaust of T-cells and CD8 cells that can lead to respiratory failure and death (Newton, Cardani, & Braciale, 2016).  Severity of the disease can lead to the hospitalization rate of 67.9 per 100,000 (CDC, 2020). COVID-19 affects individuals of age groups but particularly affects people who are 65 years and old, and immunocompromised. The overall mortality of the diseases is about 2-3%.  Virus is transmitted through respiratory droplets from expelled from a cough or sneeze and inhaled by individuals who are around the person (CDC, 2020).  The virus can also thrive on inanimate objects for 2 hours to 9 days, and can remain infections(CDC, 2020). Some common symptoms of COVID-19 are cough, shortness of breath, and fever (CDC, 2020). Individuals who are newly infected with COVID-19 will remain asymptomatic for about 2 days to 2 weeks before showing signs and symptoms, however, are still able to infect individuals around them(CDC, 2020). In San Bernardino County, there has been about more than 3,000 cases and 157 deaths (SBCPHD, 2020). 

Homelessness in San Bernardino

Since 2018, there has been a 23% increase in homelessness in San Bernardino County (PINTC, 2019)Homeless individuals had increased rates and prevalence of infectious diseases more specifically communicable acquired infectious diseases and more likely to participate in adverse behaviors such as substance abuse (alcoholism, smoking) (Badiaga, Raoult, & Brouqui, 2008). About 40% of homeless individuals report to have a chronic health conditions and 5 times more likely to be admitted to hospitalized and longer stay in the stay hospital due to these conditions compared to non-homeless individuals (Stafford & Wood, 2017). Homeless individuals have less access to health care and are less likely to test for these adverse health conditions (Stafford & Wood, 2017).   Lack of housing and usage of facilities that allows close contact (shelters, soup kitchens etc.) increases risk and prevalence for infectious disease outbreak (Stafford & Wood, 2017).  Because of these factors, there is a high percentage of homeless individuals who have compromised immune systems (Lima et al., 2020).  Homeless individuals lack of access to healthcare and proper sleeping environment might increase infection and spread of COVID-19 (Mosites, 2020). With less access to healthcare and education, people experiencing homelessness might be unaware if they are infected with COVID-19 and might unknowingly spread infection because of the lack of public safety protocols such as social distancing and wearing a mask (Tsai & Wilson, 2020). COVID-19 infectious asymptomatic individuals can easily spread through many homeless encampments and homeless individuals are much more prone to severe infection of the disease because they are much older and sicker compared to the general population. 
The purpose of this report is to determine the number of COVD-19 infected cases among the San Bernardino County unsheltered homeless population from the months of March through May 2020 using two projection models. The Loma Linda SPH team projected the worst-case scenario for the homeless population infection rate in San Bernardino County, CA given the different estimates of the unsheltered homeless population. 

# Project Description 

In coordination and updates with SBDPH CDS, projections were made using two different population estimates with assumptions about homeless health conditions for two different projection models. This study estimates the total number of infected individuals using projection models based on a “date of onset” probability trend that assumes the counts will follow a model epidemiologic infection curve. 

Population Estimates

The first estimate was the Point in Time count (PIT) from 2019 which was 1,878 for the entire county’s unsheltered homeless population ("San Bernardino Homeless Partnership 2019 San Bernardino County Homeless Count and Subpopulation Survey Final Report,"). The SB County staff considered the limitations of the PIT count and estimated the number for 2020 to be much higher at 5864, while an independent University of Pennsylvania (Penn) study estimated the number from the PIT to be 2823 (Culhane, Treglia, & Steif, 2020). We used the lowest (1,878) and highest (5,864) numbers for our projections. 

Homelessness Health Conditions

The University of Pennsylvania study estimated that homeless are in an elevated COVID-19 infection risk category. This is because the socially disenfranchised homeless experience a higher risk to COVID-19 infection due to several factors including: inadequate access to hygiene and sanitation, difficult early detection strategies in a population isolated from health care, and a high susceptibility to the infection due to their “advanced age and accelerated physical decline and mental weathering(Speer, 2016)”.

Two Models 

The two models included here are the Penn model and a stochastic model that uses the EpiModel library for R (Table 1). The Penn model uses nonlinear regression techniques to approximate rates of infection for homeless populations. This study takes those projections and calculates infection rates for the larger populations based on a direct rate conversion using their infection per population percentages (Table 1)(Culhane et al., 2020). The R EpiModel uses stochastic individual compartmental models (ICMs) on R with several input assumptions. The assumptions for this model are in Table 2. This Stochastic model expands the Susceptible, Infected and Recovered (SIR) EpiModel into more detail (Table 2).  Both models presented in this report used a “date of onset” of early March assuming no public health intervention to prevent the spread of COVID-19. They give an estimated number of infected homeless and a total number of fatalities. These models will always project an increase in cases based on an assumed “date of onset” even if the current date is beyond the peak of the curve(Churches, 2020). The EpiModel used in this study is similar to the CHIME model that assumes an exponential growth. 

Other Simulations

The other simulations ran for further investigation was Frequency Distribution(Graph 1), Baseline Simulation without safety protocols (Graph 2), Ro values (Graph 3), Baseline vs Stay at Home Simulation(Graph 4), Baseline vs Stay Home Simulation (Graph 5), Baseline vs Social Distancing at Day 15(Graph 6), Baseline vs Social Distancing at Day 30(Graph 7), Baseline vs Social Distancing(Day 15) and Stay at Home Simulation (Graph 8). Results from the simulation showed that increase public health protocols such increased social distancing at the early days of the pandemic and self-isolation would decreased the amount of individuals that are infected/ asymptomatic, case fatality, hospitalization, and decreased infected/infectious individuals. Simulations of public health interventions such as social distancing and self-isolation are seen to help flattening the curve according to the R EpiModel. As more individuals self-isolate due to social distancing and adhering to stay at home orders, the number of infectious, infected, and case fatalities decrease significantly according to simulation graphs. In our simulation table 2, starting with 5864, we have 6 number of infections in that populations by day 100 we have 4778 individuals infected (81%). For hospitalization, we have 1943 people infected (33%), 317 fatalities (.05%). 

# Results

Using public health intervention decreases the amount of infected cases and fatalities. It is significant to note, that public interventions such as social distancing and self-isolation decreases the amount of case(CDC, 2020). Many homelessness individuals do have access to adequate healthcare interventions and can spread the amount of cases COVID-19 among the nearby population in San Bernardino(Mosites, 2020). However, when comparing the simulation results and the actual cases from San Bernardino county for Homeless individuals. There were only 7 actual infected individuals in the homeless population for COVID-19(SBCPHD, 2020). This maybe due to San Bernardino County location of being more widespread in infrastructure compared to placed like New York City, where personal contact is in much more close-vicinity for the virus to spread in the area in order to infect nearby individuals(CDC, 2020). However, the information for the COVID-19 is still growing and there is much more to investigate to understand the virus. 

# Limitations 

Some of the limitations of the SEIQHFR model to assess the infection rate and movement of COVID-19 through the San Bernardino County Homeless populations in San Bernardino. Since the beginning of the pandemic, there has been a lack of testing for Covid-19 in the State of California. Because of this, it is hard to differentiate who is infected and who is not asymptomatic based on the data alone. Also, testing protocols such as drive-thru testing, only testing for COVID-19 symptoms etc. may not apply for the unsheltered homeless populations in San Bernardino and in such give inaccuracies to the actuality of initial infected cases to determine the actual infection rate.  Also, data from simulation over-estimated the number of cases compared to San Bernardino County is not match which might be dysfunction in the formal code in R or lack of testing San Bernardino County(SBCPHD, 2020).  

# Conclusions

Both models estimate the COVID-19 infection for the varying population estimates of homeless individuals in San Bernardino County. Both models produce a high estimate with the Penn model showing that 40% of all homeless are infected while the Stochastic model shows that over 81% of the homeless are infected. A study of homeless shelters in large cities showed that the COVID-19 infection percentage was as high as 66% in San Francisco with other shelters having an average of 25% infected (Speer, 2016). The high numbers reflected in these projections may truly represent a worst-case situation for homeless confined in shelters, but do not represent the true scenario of the unsheltered in SB County. However, simulation of the models overestimated and does not reflect the current situation of COVID-19 cases and deaths in San Bernardino County.  As of June 2020, there has been about 3,000 cases and 150 deaths in San Bernardino County. However, following the Public health interventions such as social distancing and self-isolation may be important for areas where COVID-19 has affected the most in the United States. Recommendations or policies to increase COVID-19 education, resources, and intervention measures for the homeless population would help address future outbreaks of in COVID-19 afflicted areas (CDC, 2020).








