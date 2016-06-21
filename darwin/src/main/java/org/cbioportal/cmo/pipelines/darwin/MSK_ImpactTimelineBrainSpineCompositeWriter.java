/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cbioportal.cmo.pipelines.darwin;

import org.cbioportal.cmo.pipelines.darwin.model.TimelineBrainSpineComposite;
//import org.cbioportal.cmo.pipelines.darwin.model.MSK_ImpactTimelineBrainSpine;
//import org.cbioportal.cmo.pipelines.darwin.MSK_ImpactTimelineBrainSpineSpecimenWriter;

import org.springframework.batch.item.*;
import org.springframework.batch.item.file.*;
import org.springframework.core.io.*;
import org.springframework.beans.factory.annotation.Value;

import java.util.*;
import org.springframework.batch.item.support.CompositeItemWriter;
/**
 *
 * @author jake
 */
public class MSK_ImpactTimelineBrainSpineCompositeWriter implements ItemStreamWriter<String>{
    
    private final FlatFileItemWriter<String> flatFileItemWriter = new FlatFileItemWriter<>();
    
    @Value("#{jobParameters[stagingDirectory]}")
    private String stagingDirectory;
    
    @Value("${darwin.timeline_bs_status}")
    private String statusFilename;
    
    @Value("${darwin.timeline_bs_specimen}")
    private String specimenFilename;
    
    @Value("${darwin.timeline_bs_treatment}")
    private String treatmentFilename;
    
    @Value("${darwin.timeline_bs_imaging}")
    private String imagingFilename;
    
    private String statusFile;
    private String specimenFile;
    private String treatmentFile;
    private String imagingFile;
    
    List<ItemStreamWriter> delegates = new ArrayList<>();
    MSK_ImpactTimelineBrainSpineStatusWriter writer1 = new MSK_ImpactTimelineBrainSpineStatusWriter();
    MSK_ImpactTimelineBrainSpineSpecimenWriter writer2 = new MSK_ImpactTimelineBrainSpineSpecimenWriter();
    MSK_ImpactTimelineBrainSpineTreatmentWriter writer3 = new MSK_ImpactTimelineBrainSpineTreatmentWriter();
    MSK_ImpactTimelineBrainSpineImagingWriter writer4 = new MSK_ImpactTimelineBrainSpineImagingWriter();
    CompositeItemWriter compWriter = new CompositeItemWriter();
    @Override
    public void close() throws ItemStreamException{
        writer1.close();
        writer2.close();
        writer3.close();
        writer4.close();
        compWriter.close();
    }
        
    @Override
    public void open(ExecutionContext executionContext) throws ItemStreamException{
        if(stagingDirectory.endsWith("/")){
            statusFile = stagingDirectory + statusFilename;
            specimenFile = stagingDirectory + specimenFilename;
            treatmentFile = stagingDirectory + treatmentFilename;
            imagingFile = stagingDirectory + imagingFilename;
        }
        else{
            statusFile = stagingDirectory + "/" + statusFilename;
            specimenFile = stagingDirectory + "/" + specimenFilename;
            treatmentFile = stagingDirectory + "/" + treatmentFilename;
            imagingFile = stagingDirectory + "/" + imagingFilename;
        }
        
        writer1.setStagingFile(statusFile);
        writer2.setStagingFile(specimenFile); 
        writer3.setStagingFile(treatmentFile);
        writer4.setStagingFile(imagingFile);
        writer1.open(executionContext);
        writer2.open(executionContext);
        writer3.open(executionContext);
        writer4.open(executionContext);
    }
    
    @Override
    public void update(ExecutionContext executionContext) throws ItemStreamException{}
    
    @Override
    public void write(List<? extends String> items) throws Exception{
        delegates.clear();
        delegates.add(writer1);
        delegates.add(writer2);
        delegates.add(writer3);
        delegates.add(writer4);
        compWriter.setDelegates(delegates);
        compWriter.write(items);
    }
    
}
