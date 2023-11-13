package org.cbioportal.cmo.pipelines.cvr.cdm;

import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.*;
import java.util.*;

import io.minio.messages.Upload;
import org.apache.log4j.Logger;
import org.cbioportal.cmo.pipelines.cvr.CVRUtilities;
import org.cbioportal.cmo.pipelines.cvr.CvrSampleListUtil;
import org.cbioportal.cmo.pipelines.cvr.clinical.CDMClinicalFieldSetMapper;
import org.cbioportal.cmo.pipelines.cvr.model.staging.CDMClinicalRecord;
import org.springframework.batch.core.StepContribution;
import org.springframework.batch.core.scope.context.ChunkContext;
import org.springframework.batch.core.step.tasklet.Tasklet;
import org.springframework.batch.item.ItemStreamException;
import org.springframework.batch.item.file.FlatFileHeaderCallback;
import org.springframework.batch.item.file.FlatFileItemReader;
import org.springframework.batch.item.file.FlatFileItemWriter;
import org.springframework.batch.item.file.mapping.DefaultLineMapper;
import org.springframework.batch.item.file.transform.DelimitedLineTokenizer;
import org.springframework.batch.item.file.transform.PassThroughLineAggregator;
import org.springframework.batch.repeat.RepeatStatus;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.core.io.FileSystemResource;
import io.minio.*;

/**
 *
 * @author averyniceday
 *
 * CDMPublisherTasklet reads in clinical file and constructs a new file that contains the CDM information
 */
public class CDMPublisherTasklet implements Tasklet {
    @Value("#{jobParameters[testingMode]}")
    private boolean testingMode;

    @Value("#{jobParameters[stagingDirectory]}")
    private String stagingDirectory;

    @Value("#{jobParameters[clinicalFilename]}")
    private String clinicalFilename;

    @Autowired
    private CvrSampleListUtil cvrSampleListUtil;

    @Autowired
    private CVRUtilities cvrUtilities;

    @Autowired
    private MinioClient cdmMinioClient;

    private final ObjectMapper mapper = new ObjectMapper();
    private final Logger log = Logger.getLogger(CDMPublisherTasklet.class);

    @Override
    public RepeatStatus execute(StepContribution sc, ChunkContext cc) throws Exception {
        if (testingMode) {
            log.info("[TEST MODE] samples will not be published to smile");
            //return RepeatStatus.FINISHED;
        }
        // iterate through json file and if sample id is in list of samples that were
        // consumed successfully then publish its metadata to smile
        File CDMFile = File.createTempFile("cdm", ".txt");
        String CDMFilepath = CDMFile.getAbsolutePath();
        processCDMClinicalRecords(sc, clinicalFilename, CDMFilepath);

        log.info("temp file located at " + CDMFilepath);
        //connect and publish file to minio
        
        System.out.println("UPLOADING " + CDMFilepath);
        ObjectWriteResponse resp = cdmMinioClient.uploadObject(UploadObjectArgs.builder()
                .bucket("cdm-data")
                .object("cbioportal/CDMData.txt")
                .filename(CDMFilepath)
                .contentType("text/plain")
                .build());
        return RepeatStatus.FINISHED;
    }

    private void processCDMClinicalRecords(StepContribution sc, String clinicalFilename, String CDMFilepath) throws Exception {
        File mskimpactClinicalFile = new File(stagingDirectory, clinicalFilename);
        if (!mskimpactClinicalFile.exists()) {
            log.error("File does not exist - skipping data loading from clinical file: " + mskimpactClinicalFile.getName());
            return;
        }
        log.info("Loading clinical data from: " + mskimpactClinicalFile.getName());
        DelimitedLineTokenizer tokenizer = new DelimitedLineTokenizer(DelimitedLineTokenizer.DELIMITER_TAB);
        DefaultLineMapper<CDMClinicalRecord> mapper = new DefaultLineMapper<>();
        mapper.setLineTokenizer(tokenizer);
        mapper.setFieldSetMapper(new CDMClinicalFieldSetMapper());

        FlatFileItemReader<CDMClinicalRecord> reader = new FlatFileItemReader<>();
        reader.setResource(new FileSystemResource(mskimpactClinicalFile));
        reader.setLineMapper(mapper);
        reader.setLinesToSkip(1);
        reader.open(sc.getStepExecution().getJobExecution().getExecutionContext());

        List<String> cdmRecords = new ArrayList<String>();
        CDMClinicalRecord record;
        try {
            while ((record = reader.read()) != null) {
                List<String> to_add = new ArrayList<String>();
                for (String field : CDMClinicalRecord.getFieldNames()) {
                    try {
                        to_add.add(cvrUtilities.convertWhitespace(record.getClass().getMethod("get" + field).invoke(record).toString().trim()));
                    } catch (NullPointerException e) {
                        log.error("Null pointer expection: " + field);
                        to_add.add("");
                    }
                } 
                cdmRecords.add(String.join("\t", to_add));
            }
        }
        catch (Exception e) {
            log.error("Error reading data from clinical file: " + mskimpactClinicalFile.getName());
            throw new ItemStreamException(e);
        }
        finally {
            reader.close();
        }

        FlatFileItemWriter<String> writer = new FlatFileItemWriter<>();
        PassThroughLineAggregator aggr = new PassThroughLineAggregator();
        writer.setResource(new FileSystemResource(CDMFilepath));
        writer.setLineAggregator(aggr);
        writer.setTransactional(false);
        writer.setForceSync(true);
        writer.setHeaderCallback(new FlatFileHeaderCallback() {
            @Override
            public void writeHeader(Writer writer) throws IOException {
                writer.write(String.join("\t", CDMClinicalRecord.getFieldNames()));
            }
        });
        writer.setResource(new FileSystemResource(CDMFilepath));
        writer.open(sc.getStepExecution().getJobExecution().getExecutionContext());
        try {
            writer.write(cdmRecords);
            writer.close();
        } catch (Exception e) {
            log.error("Error writing data to clinical file: " + CDMFilepath);
            throw new ItemStreamException(e);
        }
        finally {
            writer.close();
            System.out.println("WTF");
        }
    }

}
