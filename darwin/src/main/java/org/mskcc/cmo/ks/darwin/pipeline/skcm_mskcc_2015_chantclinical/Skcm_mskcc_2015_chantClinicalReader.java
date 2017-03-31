/*
 * Copyright (c) 2016 Memorial Sloan-Kettering Cancer Center.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF MERCHANTABILITY OR FITNESS
 * FOR A PARTICULAR PURPOSE. The software and documentation provided hereunder
 * is on an "as is" basis, and Memorial Sloan-Kettering Cancer Center has no
 * obligations to provide maintenance, support, updates, enhancements or
 * modifications. In no event shall Memorial Sloan-Kettering Cancer Center be
 * liable to any party for direct, indirect, special, incidental or
 * consequential damages, including lost profits, arising out of the use of this
 * software and its documentation, even if Memorial Sloan-Kettering Cancer
 * Center has been advised of the possibility of such damage.
 */

/*
 * This file is part of cBioPortal CMO-Pipelines.
 *
 * cBioPortal is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
package org.mskcc.cmo.ks.darwin.pipeline.skcm_mskcc_2015_chantclinical;

import org.mskcc.cmo.ks.darwin.pipeline.mskimpactdemographics.MskimpactPatientDemographicsReader;
import org.mskcc.cmo.ks.darwin.pipeline.model.Skcm_mskcc_2015_chantClinicalRecord;
import java.util.*;
import org.springframework.batch.item.*;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.transaction.annotation.Transactional;
import org.apache.log4j.Logger;
import com.querydsl.core.types.Projections;
import static com.querydsl.core.alias.Alias.*;
import static com.querydsl.core.alias.Alias.alias;
import com.querydsl.sql.SQLQueryFactory;
import java.lang.reflect.Method;
import org.apache.commons.lang.StringUtils;

/**
 *
 * @author heinsz
 */
public class Skcm_mskcc_2015_chantClinicalReader implements ItemStreamReader<Skcm_mskcc_2015_chantClinicalRecord>{
    
    @Value("${darwin.skcm_mskcc_2015_chant.staging_path_current_view}")
    private String stagingPathCurrentView;
    
    @Value("${darwin.skcm_mskcc_2015_chant.primary_view}")
    private String primaryView;

    @Value("${darwin.skcm_mskcc_2015_chant.general_view}")
    private String generalView;     
    
    @Value("${darwin.skcm_mskcc_2015_chant.impact_view}")
    private String impactView;
    
    @Value("${darwin.skcm_mskcc_2015_chant.metastat_sites_view}")
    private String metastatSitesView;      
    
    @Autowired
    SQLQueryFactory darwinQueryFactory;
    
    Logger log = Logger.getLogger(MskimpactPatientDemographicsReader.class);
    
    List<Skcm_mskcc_2015_chantClinicalRecord> melanomaClinicalRecords;
    
    @Override
    public void open(ExecutionContext executionContext) throws ItemStreamException {
        this.melanomaClinicalRecords = getMelanomaClinicalRecords();
    }
    
    @Override
    public void update(ExecutionContext executionContext) throws ItemStreamException{}
    
    @Override
    public void close() throws ItemStreamException{}
    
    @Override
    public Skcm_mskcc_2015_chantClinicalRecord read() throws Exception{
        if(!melanomaClinicalRecords.isEmpty()) {
            return melanomaClinicalRecords.remove(0);
        }
        return null;
    }
    
    private List<Skcm_mskcc_2015_chantClinicalRecord> getMelanomaClinicalRecords() {
        log.info("Start of Darwin skcm_mskcc_2015_chant clinical records query..");
        List<Skcm_mskcc_2015_chantClinicalRecord> rawDarwinMelanomaRecords = getRecordsFromDatabase();
        List<Skcm_mskcc_2015_chantClinicalRecord> combinedRecords = combineRecords(rawDarwinMelanomaRecords);
        return combinedRecords;
    }
    
    @Transactional
    private List<Skcm_mskcc_2015_chantClinicalRecord> getRecordsFromDatabase() {
        Skcm_mskcc_2015_chantClinicalRecord qStagingPathCurrentView = alias(Skcm_mskcc_2015_chantClinicalRecord.class, stagingPathCurrentView);
        Skcm_mskcc_2015_chantClinicalRecord qPrimaryView = alias(Skcm_mskcc_2015_chantClinicalRecord.class, primaryView);      
        Skcm_mskcc_2015_chantClinicalRecord qGeneralView = alias(Skcm_mskcc_2015_chantClinicalRecord.class, generalView);
        Skcm_mskcc_2015_chantClinicalRecord qImpactView = alias(Skcm_mskcc_2015_chantClinicalRecord.class, impactView);
        Skcm_mskcc_2015_chantClinicalRecord qMetastatSitesView = alias(Skcm_mskcc_2015_chantClinicalRecord.class, metastatSitesView);        
        
        List<Skcm_mskcc_2015_chantClinicalRecord> darwinMelanomaRecords = darwinQueryFactory.selectDistinct(Projections.constructor(Skcm_mskcc_2015_chantClinicalRecord.class,
                $(qStagingPathCurrentView.getMELSPC_PTID()),
                $(qStagingPathCurrentView.getMELSPC_STAGE_YEAR()),
                $(qStagingPathCurrentView.getMELSPC_STG_GRP_NAME()),
                $(qGeneralView.getMELG_PTID()),
                $(qGeneralView.getMELG_STS_DESC()),
                $(qGeneralView.getMELG_STS_SRCE_DESC()),
                $(qGeneralView.getMELG_ACTV_STS_DESC()),
                $(qGeneralView.getMELG_DERMAGRPHX_DESC()),
                $(qGeneralView.getMELG_PRES_STG_YEAR()),
                $(qGeneralView.getMELG_FAMILY_HX_DESC()),
                $(qGeneralView.getMELG_1ST_RECUR_YEAR()),
                $(qGeneralView.getMELG_LOCAL_DESC()),
                $(qGeneralView.getMELG_NODAL_DESC()),
                $(qGeneralView.getMELG_INTRANSIT_DESC()),
                $(qGeneralView.getMELG_SYS_DESC()),
                $(qGeneralView.getMELG_RECUR_NDSZ_DESC()),
                $(qGeneralView.getMELG_RECUR_NODAL_NO()),
                $(qGeneralView.getMELG_LDH()),
                $(qGeneralView.getMELG_LDH_YEAR()),
                $(qGeneralView.getMELG_METS_DESC()),
                $(qGeneralView.getMELG_ADJVNT_TX_DESC()),
                $(qGeneralView.getMELG_SYS_TX_DESC()),
                $(qGeneralView.getMELG_RAD_TX_DESC()),
                $(qGeneralView.getMELG_SURG_DESC()),
                $(qGeneralView.getMELG_TISSUE_BANK_AVAIL()),
                $(qPrimaryView.getMELP_PTID()),
                $(qPrimaryView.getMELP_PRIM_SEQ()),
                $(qPrimaryView.getMELP_DX_YEAR()),
                $(qPrimaryView.getMELP_MSK_REVIEW_DESC()),
                $(qPrimaryView.getMELP_THICKNESS_MM()),
                $(qPrimaryView.getMELP_CLARK_LVL_DESC()),
                $(qPrimaryView.getMELP_ULCERATION_DESC()),
                $(qPrimaryView.getMELP_SITE_DESC()),
                $(qPrimaryView.getMELP_SUB_SITE_DESC()),
                $(qPrimaryView.getMELP_TILS_DESC()),
                $(qPrimaryView.getMELP_REGRESSION_DESC()),
                $(qPrimaryView.getMELP_MARGINS_DESC()),
                $(qPrimaryView.getMELP_MITIDX_UNK_DESC()),
                $(qPrimaryView.getMELP_HIST_TYPE_DESC()),
                $(qPrimaryView.getMELP_SATELLITES_DESC()),
                $(qPrimaryView.getMELP_EXT_SLIDES_DESC()),
                $(qPrimaryView.getMELP_LNORG_DX_DESC()),
                $(qPrimaryView.getMELP_LNCLIN_STS_DESC()),
                $(qPrimaryView.getMELP_LNSENTINBX_DESC()),
                $(qPrimaryView.getMELP_LNSENTINBX_YEAR()),
                $(qPrimaryView.getMELP_LNPROLYSCT_DESC()),
                $(qPrimaryView.getMELP_LNPROSUCC_DESC()),
                $(qPrimaryView.getMELP_LNDSCT_CMP_DESC()),
                $(qPrimaryView.getMELP_LNDSCT_YEAR()),
                $(qPrimaryView.getMELP_LNMATTED_DESC()),
                $(qPrimaryView.getMELP_LNEXTNODST_DESC()),
                $(qPrimaryView.getMELP_LNINTRMETS_DESC()),
                $(qPrimaryView.getMELP_LNSIZE()),
                $(qPrimaryView.getMELP_LNSIZE_UNK_DESC()),
                $(qPrimaryView.getMELP_LNSLNLARG_SIZE()),
                $(qPrimaryView.getMELP_LNIHC_DESC()),
                $(qPrimaryView.getMELP_LNIMM_S100_DESC()),
                $(qPrimaryView.getMELP_LNIMMHMB45_DESC()),
                $(qPrimaryView.getMELP_LNIMM_MELA_DESC()),
                $(qImpactView.getMELI_PTID()),
                $(qImpactView.getMELI_DMP_PATIENT_ID()),
                $(qImpactView.getMELI_DMP_SAMPLE_ID()),
                $(qImpactView.getMELI_REPORT_YEAR()),
                $(qImpactView.getMELI_PROCEDURE_YEAR()),
                $(qImpactView.getMELI_TUMOR_TYPE()),
                $(qImpactView.getMELI_PRIMARY_SITE()),
                $(qImpactView.getMELI_MET_SITE()),
                $(qMetastatSitesView.getMELMS_PTID()),
                $(qMetastatSitesView.getMELMS_SITE_TYPE_DESC()),
                $(qMetastatSitesView.getMELMS_SITE_DESC()),
                $(qMetastatSitesView.getMELMS_SITE_YEAR())))
                .from($(qStagingPathCurrentView))
                .fullJoin($(qPrimaryView))
                .on($(qStagingPathCurrentView.getMELSPC_PTID()).eq($(qPrimaryView.getMELP_PTID())))
                .fullJoin($(qGeneralView))
                .on($(qStagingPathCurrentView.getMELSPC_PTID()).eq($(qGeneralView.getMELG_PTID())).or($(qPrimaryView.getMELP_PTID()).eq($(qGeneralView.getMELG_PTID()))))
                .fullJoin($(qImpactView))
                .on($(qStagingPathCurrentView.getMELSPC_PTID()).eq($(qImpactView.getMELI_PTID())).or($(qPrimaryView.getMELP_PTID()).eq($(qImpactView.getMELI_PTID())).or($(qGeneralView.getMELG_PTID()).eq($(qImpactView.getMELI_PTID())))))
                .fullJoin($(qMetastatSitesView))
                .on($(qGeneralView.getMELG_PTID()).eq($(qMetastatSitesView.getMELMS_PTID())).or($(qPrimaryView.getMELP_PTID()).eq($(qMetastatSitesView.getMELMS_PTID()))).or($(qStagingPathCurrentView.getMELSPC_PTID()).eq($(qMetastatSitesView.getMELMS_PTID()))))
                .fetch();
        return darwinMelanomaRecords;
    }
    
    private List<Skcm_mskcc_2015_chantClinicalRecord> combineRecords(List<Skcm_mskcc_2015_chantClinicalRecord> records) {
        Map<String, Skcm_mskcc_2015_chantClinicalRecord> recordMap = new HashMap();
        List<Skcm_mskcc_2015_chantClinicalRecord> combinedRecords = new ArrayList<>();
        
        for (Skcm_mskcc_2015_chantClinicalRecord record : records) {
            if (!recordMap.containsKey(record.getSAMPLE_ID())) {
                recordMap.put(record.getSAMPLE_ID(), record);
            }
            else {
                Skcm_mskcc_2015_chantClinicalRecord existingRecord = recordMap.get(record.getSAMPLE_ID());
                recordMap.put(record.getSAMPLE_ID(), combineRecords(existingRecord, record));
            }
        }
        combinedRecords.addAll(recordMap.values());
        return combinedRecords;
    }
    
    private Skcm_mskcc_2015_chantClinicalRecord combineRecords(Skcm_mskcc_2015_chantClinicalRecord record1, Skcm_mskcc_2015_chantClinicalRecord record2) {
        Skcm_mskcc_2015_chantClinicalRecord combined = new Skcm_mskcc_2015_chantClinicalRecord();
        try {
            for (String field : combined.getAllVariables()) {
                    Method fieldGetter = combined.getClass().getMethod("get" + field);
                    List<String> values = new ArrayList();
                    values.add((String) fieldGetter.invoke(record1));
                    values.add((String) fieldGetter.invoke(record2));
                    combined.getClass().getMethod("set" + field, String.class).invoke(combined, StringUtils.join(values, "|"));
            }
        }
        catch (Exception e) {
            log.error("Failed to combine records!" + e.getMessage());
            throw new ItemStreamException(e);        
        }
        return combined;
    }
}