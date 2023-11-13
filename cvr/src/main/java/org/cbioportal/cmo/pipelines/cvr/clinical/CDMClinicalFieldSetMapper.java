package org.cbioportal.cmo.pipelines.cvr.clinical;

import org.apache.log4j.Logger;
import org.cbioportal.cmo.pipelines.cvr.model.staging.CDMClinicalRecord;
import org.springframework.batch.item.file.mapping.FieldSetMapper;
import org.springframework.batch.item.file.transform.FieldSet;
import org.springframework.validation.BindException;

import java.util.List;

public class CDMClinicalFieldSetMapper implements FieldSetMapper<CDMClinicalRecord> {
    Logger log = Logger.getLogger(CVRClinicalFieldSetMapper.class);

    @Override
    public CDMClinicalRecord mapFieldSet(FieldSet fs) throws BindException {
        CDMClinicalRecord record = new CDMClinicalRecord();
        List<String> fields = CDMClinicalRecord.getFieldNames();

        for (int i = 0; i < fields.size(); i++) {
            String field = fields.get(i);
            try {
                record.getClass().getMethod("set" + field, String.class).invoke(record, fs.readString(i));
            } catch (Exception e) {
                if (e.getClass().equals(NoSuchMethodException.class)) {
                    log.info("No set method exists for " + field);
                }
            }
        }
        return record;
    }
}
