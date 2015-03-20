// -------------------------------------------------------------------
// Biomedical Engineering Program
// Department of Systems Design Engineering
// University of Waterloo
//
// Student Name:     <JOYCE ZHANG>
// Userid:           <J495ZHAN>
//
// Assignment:       <PA2>
// Submission Date:  <October 6th, 2014>
// 
// I declare that, other than the acknowledgements listed below, 
// this program is my original work.
//
// Acknowledgements:
// <N/A>
// -------------------------------------------------------------------

using System;
static class geneticDriftSimulator
{
	static void Main( )
	{
        // Generates random numbers for number of people with gene of interest.
        Random r = new Random( );
       
        //Collecting input: 
        Console.Write("Enter the number of experiments: ");
        int nExperiments = int.Parse(Console.ReadLine());

        
        Console.Write("Enter the population size: ");
        int peoplePerGen = int.Parse(Console.ReadLine());

        
        Console.Write("Enter the initial number holding the gene: ");
        int initialWithGene = int.Parse(Console.ReadLine());
        
        //Declaring variables
        int ExperimentNumber = 0; // the experiment number being tested
        int noGeneStableState = 0; // a stable state without any of the gene
        int allGeneStableState = 0; // a stable state with all with the gene
        
        int withGene = 0; // an individual with the gene
        int withoutGene = 0; // an individual without the gene
        
        // number of people in population with the gene
        int withGenePop = initialWithGene; 
        
        int totalGen = 0; // total generations tested in all the experiments
        
   
        //This loop runs the correct number of experiments
        while (ExperimentNumber < nExperiments )
        {    
            
            //resets the variables for each experiment
            withGenePop = initialWithGene;
            withoutGene = peoplePerGen - initialWithGene;
        
            //This loop tests when the generation is stable.
            while( (withGenePop < peoplePerGen) && (withoutGene < peoplePerGen))
                {
                    // resetting the variables
                    withGene = 0;
                    withoutGene = 0;
                    
                    //counts the individuals which have been tested for gene in 
                    //the same generation
                    int indivCounter = 0;
                    
                    //This loop tests each individual for the gene
                    while (indivCounter < peoplePerGen)
                    {
                        if ( r.Next(peoplePerGen) < withGenePop) 
                        {
                            withGene++;
                  
                        }
                        else
                        {
                            withoutGene++;
                      
                        }
                         indivCounter++;
                    }
                    
                    // Updates for number people in the current generation with
                    //the gene.
                    withGenePop = withGene;
                    
                    // collects the total number of generations for avgGen
                    totalGen++; 
                }
            
            //This records the number of each stable state
            //reached in the experiments
            if (withGenePop == peoplePerGen)
            {
                allGeneStableState++;
            }
            else if (withoutGene == peoplePerGen)
            {
                noGeneStableState++;
            }
            
            ExperimentNumber++;
        }
        
        //Conversions and calculations for avgGen:
        float floatTotalGen = (float) totalGen;
        float floatNExperiments = (float) nExperiments;
        float avgGen = floatTotalGen / floatNExperiments;

        // Final output to user:
        Console.WriteLine( "Experiments ending with none holding the gene: {0} "
            , allGeneStableState);
        
        Console.WriteLine( "Experiments ending this all holding the gene: {0} "
            , noGeneStableState);
        
        Console.WriteLine( "Average generations to reach a stable state: {0:F3}"
            , avgGen );  
    }
}    
        
        