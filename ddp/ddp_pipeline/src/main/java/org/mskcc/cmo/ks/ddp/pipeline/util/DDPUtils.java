/*
 * Copyright (c) 2018 Memorial Sloan-Kettering Cancer Center.
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

package org.mskcc.cmo.ks.ddp.pipeline.util;

import org.mskcc.cmo.ks.ddp.source.composite.DDPCompositeRecord;
import org.mskcc.cmo.ks.ddp.source.model.PatientDiagnosis;

import com.google.common.base.Strings;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;
import org.apache.commons.lang.StringUtils;

/**
 *
 * @author ochoaa
 */
public class DDPUtils {
    public static final Double DAYS_TO_YEARS_CONVERSION = 365.2422;
    public static final Double DAYS_TO_MONTHS_CONVERSION = 30.4167;
    public static final List<String> NULL_EMPTY_VALUES = Arrays.asList(new String[]{"NA", "N/A"});

    /**
     * Resolve and anonymize patient current age.
     *
     * If birth year is null/empty then return anonymized patient age.
     * Otherwise calculate from difference between reference date and birth date,
     * where reference date is either date of death if pt is deceased or current date.
     *
     * @param compositeRecord
     * @return
     * @throws ParseException
     */
    public static String resolvePatientCurrentAge(DDPCompositeRecord compositeRecord) throws ParseException {
        // if patient birth date is null/empty then use current age value
        if (Strings.isNullOrEmpty(compositeRecord.getPatientBirthDate())) {
            return anonymizePatientAge(compositeRecord.getPatientAge());
        }
        Long birthDateInDays = getDateInDays(compositeRecord.getPatientBirthDate());
        Long referenceDateInDays;
        // use current date as reference date if patient not deceased, otherwise use date of death
        if (!Strings.isNullOrEmpty(compositeRecord.getPatientDeathDate())) {
            referenceDateInDays = getDateInDays(compositeRecord.getPatientDeathDate());
        }
        else {
            referenceDateInDays = getDateInDays(new Date());
        }
        Double age = (referenceDateInDays - birthDateInDays) / (DAYS_TO_YEARS_CONVERSION);
        return anonymizePatientAge(age.intValue());
    }

    /**
     * Resolve and anonymize patient age at diagnosis.
     *
     * If birth year is null/empty then return anonymized patient age.
     * Otherwise calculate from difference between current date and birth date.
     *
     * @param compositeRecord
     * @return
     * @throws ParseException
     */
    public static String resolvePatientAgeAtDiagnosis(DDPCompositeRecord compositeRecord) throws ParseException {
        // if patient birth date is null/empty then use current age value
        if (Strings.isNullOrEmpty(compositeRecord.getPatientBirthDate())) {
            return anonymizePatientAge(compositeRecord.getPatientAge());
        }
        Long birthDateInDays = getDateInDays(compositeRecord.getPatientBirthDate());
        Long currentDateInDays = getDateInDays(new Date());
        Double age = (currentDateInDays - birthDateInDays) / (DAYS_TO_YEARS_CONVERSION);
        return anonymizePatientAge(age.intValue());
    }

    /**
     * Returns anonymized patient age as string.
     * @param age
     * @return
     */
    private static String anonymizePatientAge(Integer age) {
        if (age >= 90) {
            age = 90;
        }
        else if (age <= 18) {
            age = 18;
        }
        return String.valueOf(age);
    }

    /**
     * Standardize patient sex value.
     *
     * @param compositeRecord
     * @return
     */
    public static String resolvePatientSex(DDPCompositeRecord compositeRecord) {
        String sex = "NA";
        if (!isNullEmptyValue(compositeRecord.getPatientSex())) {
            if (compositeRecord.getPatientSex().equalsIgnoreCase("M") ||
                    compositeRecord.getPatientSex().equalsIgnoreCase("MALE")) {
                sex = "Male";
            }
            else if (compositeRecord.getPatientSex().equalsIgnoreCase("F") ||
                    compositeRecord.getPatientSex().equalsIgnoreCase("FEMALE")){
                sex = "Female";
            }
        }
        return sex;
    }

    /**
     * Standardize OS_STATUS value.
     * If CohortPatient is null then resolve from PatientDemographics,
     * which should never be null
     *
     * @param compositeRecord
     * @return
     */
    public static String resolveOsStatus(DDPCompositeRecord compositeRecord) {
        String osStatus = "NA";
        if (compositeRecord.getCohortPatient() != null &&
                !Strings.isNullOrEmpty(compositeRecord.getCohortPatient().getPTVITALSTATUS())) {
            if (compositeRecord.getCohortPatient().getPTVITALSTATUS().equalsIgnoreCase("ALIVE")) {
                osStatus = "LIVING";
            }
            else {
                osStatus = "DECEASED";
            }
        }
        else if (compositeRecord.getPatientDemographics() != null) {
            if (!Strings.isNullOrEmpty(compositeRecord.getPatientDemographics().getDeceasedDate()) ||
                !Strings.isNullOrEmpty(compositeRecord.getPatientDemographics().getPTDEATHDTE())) {
                osStatus = "DECEASED";
            }
            else {
                osStatus = "LIVING";
            }
        }
        return osStatus;
    }

    /**
     * Calculate OS_MONTHS.
     *
     * Note: In some cases, patients may not have any tumor diagnoses in the system yet. These are NA.
     *
     * @param osStatus
     * @param compositeRecord
     * @return
     * @throws ParseException
     */
    public static String resolveOsMonths(String osStatus, DDPCompositeRecord compositeRecord) throws ParseException {
        String osMonths = "NA";
        Long referenceAgeInDays = (osStatus.equals("LIVING")) ?
                getDateInDays(compositeRecord.getPatientDemographics().getPLALASTACTVDTE()) :
                getDateInDays(compositeRecord.getPatientDemographics().getDeceasedDate());
        Long firstTumorDiagnosisDateInDays = getFirstTumorDiagnosisDateInDays(compositeRecord.getPatientDiagnosis());
        if (referenceAgeInDays != null && firstTumorDiagnosisDateInDays != null) {
            osMonths = String.format("%.3f", (referenceAgeInDays - firstTumorDiagnosisDateInDays) / DAYS_TO_MONTHS_CONVERSION);
        }
        return osMonths;
    }

    /**
     * Find and return the earliest patient tumor diagnosis date in days.
     *
     * @param patientDiagnosis
     * @return
     * @throws ParseException
     */
    private static Long getFirstTumorDiagnosisDateInDays(List<PatientDiagnosis> patientDiagnosis) throws ParseException {
        return getDateInDays(getFirstTumorDiagnosisDate(patientDiagnosis));
    }

    /**
     * Find and return the earliest patient tumor diagnosis date as a string.
     *
     * @param patientDiagnosis
     * @return
     * @throws ParseException
     */
    public static String getFirstTumorDiagnosisDate(List<PatientDiagnosis> patientDiagnosis) throws ParseException {
        if (patientDiagnosis == null || patientDiagnosis.isEmpty()) {
            return null;
        }
        SimpleDateFormat sfd = new SimpleDateFormat("yyyy-MM-dd");
        String firstTumorDiagnosisDateValue = null;
        Date firstTumorDiagnosisDate = null;
        for (PatientDiagnosis diagnosis : patientDiagnosis) {
            if (Strings.isNullOrEmpty(diagnosis.getTumorDiagnosisDate())) {
                continue;
            }
            if (firstTumorDiagnosisDate == null) {
                firstTumorDiagnosisDate = sfd.parse(diagnosis.getTumorDiagnosisDate());
                firstTumorDiagnosisDateValue = diagnosis.getTumorDiagnosisDate();
            }
            else {
                Date currentTumorDiagnosisDate = sfd.parse(diagnosis.getTumorDiagnosisDate());
                if (currentTumorDiagnosisDate.before(firstTumorDiagnosisDate)) {
                    firstTumorDiagnosisDate = currentTumorDiagnosisDate;
                    firstTumorDiagnosisDateValue = diagnosis.getTumorDiagnosisDate();
                }
            }
        }
        return firstTumorDiagnosisDateValue;
    }

    /**
     * Resolves the date offset in days between the timeline event and reference date (date of birth or first tumor diagnosis date).
     * @param referenceDate
     * @param eventDate
     * @return
     * @throws ParseException
     */
    public static String resolveTimelineEventDateInDays(String referenceDate, String eventDate) throws ParseException {
        Long birthDateInDays = getDateInDays(referenceDate);
        Long eventDateInDays = getDateInDays(eventDate);
        Long timelineEventInDays = (eventDateInDays - birthDateInDays);
        return timelineEventInDays.toString();
    }

    /**
     * Calculates the date in days.
     *
     * @param date
     * @return
     * @throws ParseException
     */
    private static Long getDateInDays(Date date) throws ParseException {
        if (date != null) {
            return date.getTime() / (1000 * 60 * 60 * 24); // milliseconds -> minutes -> hours -> days
        }
        return null;
    }

    /**
     * Calculates the date in days.
     *
     * @param dateValue
     * @return
     * @throws ParseException
     */
    private static Long getDateInDays(String dateValue) throws ParseException {
        if (!Strings.isNullOrEmpty(dateValue)) {
            SimpleDateFormat sfd = new SimpleDateFormat("yyyy-MM-dd");
            Date date = sfd.parse(dateValue);
            return date.getTime() / (1000 * 60 * 60 * 24); // milliseconds -> minutes -> hours -> days
        }
        return null;
    }

    /**
     * Converts record as tab-delimited string of values.
     *
     * @param object
     * @return
     * @throws Exception
     */
    @SuppressWarnings("unchecked")
    public static String constructRecord(Object object) throws Exception {
        List<String> fields = (List<String>) object.getClass().getMethod("getFieldNames").invoke(object); // unchecked cast
        List<String> record = new ArrayList<>();
        for (String field : fields) {
            String value = object.getClass().getMethod("get" + field).invoke(object).toString();
            record.add(value.trim());
        }
        return StringUtils.join(record, "\t");
    }

    /**
     * Returns whether value is null/empty string.
     * @param value
     * @return
     */
    public static Boolean isNullEmptyValue(String value) {
        return (Strings.isNullOrEmpty(value) || NULL_EMPTY_VALUES.contains(value.trim()));
    }
}