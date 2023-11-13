package org.cbioportal.cmo.pipelines.cvr.model.staging;

import java.util.ArrayList;
import java.util.List;

public class CDMClinicalRecord {
    private String PATIENT_ID;
    private String SAMPLE_ID;

    public CDMClinicalRecord() {
    }

    public String getPATIENT_ID() {
        return PATIENT_ID;
    }

    public void setPATIENT_ID(String PATIENT_ID) {
        this.PATIENT_ID = PATIENT_ID;
    }

    public String getSAMPLE_ID() {
        return SAMPLE_ID;
    }

    public void setSAMPLE_ID(String SAMPLE_ID) {
        this.SAMPLE_ID = SAMPLE_ID;
    }

    public static List<String> getFieldNames() {
        List<String> fieldNames = new ArrayList<>();
        fieldNames.add("PATIENT_ID");
        fieldNames.add("SAMPLE_ID");
        return fieldNames;
    }
}
