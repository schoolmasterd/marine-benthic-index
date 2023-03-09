OBTAIN ENVIRONMENTAL DATA

1.  **From tbl_Chemistry:** In query, specify

-   desired Project,

-   desired SamplingYear,

-   LabRep=\"1\",

-   ParameterCode=\"Total Fines\" OR \"Gravel\" OR \"Inorganic Carbon\",

-   SampleType=\"Result\",

-   QACode=\"A\".

> If not using all stations in the desired project and year, specify
> desired StationID list. If not using all field replicates, specify
> FieldReplicate=\"1\".

a)  If FieldReplicate 1 is not available or is unusable, use field split
    > if available (usually, but not always, designated as
    > FieldReplicate 4 for Long-term or FieldReplicate 2 for Urban
    > Bays).

b)  If LabRep 1 is not available or is unusable, use another lab
    > duplicate (usually LabRep 2) if available.

c)  Ensure that there is exactly one record per parameter per sample.

d)  Concatenate Project, Year, Station, and FieldReplicate, separated by
    > underscores, to create the Sample ID. Remove the space in \"Urban
    > Bays\".

e)  Rename \"Total Fines\" as \"Fines\" and \"Inorganic Carbon\" as
    > \"InorgC\". Case is important.

f)  If any Gravel is nondetect (Qualifier = U or UJ), ensure that the
    > Result value is zero.

g)  If any InorgC Result value is \<0, set it to zero.

h)  Create one sample record per Sample containing the results for each
    > parameter in separate fields.

i)  The fields required are Sample, Fines, Gravel, and InorgC.

```{=html}
<!-- -->
```
2.  **From tbl_GrabEvent:** In query, specify

-   desired Project,

-   desired SamplingYear,

-   BenthicInfauna=TRUE,

-   GrabFailCode=\"None\".

> If not using all stations in the desired project and year, specify
> desired StationID list. If not using all field replicates, specify
> FieldReplicate=\"1\".

a)  If there is not exactly one record for each sample (either too few
    > or too many), check the written field and navigation logs to
    > determine which grab\'s record to use for that sample. It may be
    > that there is an error in the database identifying which grab was
    > used for benthos or that the grab failure code is something other
    > than None or that a second benthos grab was obtained for other
    > purposes, such as DNA analysis.

b)  Concatenate Project, Year, Station, and FieldReplicate, separated by
    > underscores, to create the Sample ID. Remove the space in \"Urban
    > Bays\".

c)  If the MeterWheelDepth is blank or set to a nonsense value (e.g.,
    > -9999), obtain the StationDepth from tbl_StationOccupation to use
    > instead.

d)  Rename \"MeterWheelDepth\" as \"Depth\".

e)  If the depth values are text, convert the values to numeric.

f)  The fields required are Sample and Depth.

```{=html}
<!-- -->
```
3.  **From tbl_field notes:** In query, specify

-   desired Project,

-   desired SamplingYear.

> If not using all stations in the desired project and year, specify
> desired StationID list.

a)  If there is not exactly one record for each sample, check the
    > written field and navigation logs to determine where the error
    > lies and make corrections accordingly.

b)  If the salinity, temperature, or penetration values are text,
    > convert the values to numeric.

c)  If the salinity, temperature, or penetration values are missing or
    > nonsense values (e.g., -9999), check the written field logs to
    > determine the appropriate value (may be recorded for different
    > grab, field split, or field replicate). If the values are still
    > missing, try estimating the values as follows:

d)  Penetration: Use the penetration depth from the same station in
    > another year.

e)  Salinity or temperature: If another station nearby (same embayment)
    > was sampled within an hour or two, use the value recorded for that
    > station. For deep stations, it may be reasonable to use the
    > salinity from the same station in another year.

f)  Rename \"Sediment Temperature\" as \"Temperature\". Case is
    > important.

g)  Concatenate Project, Year, Station, and FieldReplicate, separated by
    > underscores, to create the Sample ID. Remove the space in \"Urban
    > Bays\".

h)  The fields required are Sample, Penetration, Salinity, and
    > Temperature.

```{=html}
<!-- -->
```
4.  **Combine:** Match on the Sample IDs to concatenate the records to
    create a single record for each sample with the following fields:
    Sample, Depth, Penetration, Salinity, Temperature, Fines, Gravel,
    and InorgC.

5.  Save as a csv file.

OBTAIN BENTHOS DATA

6.  **Join tbl_InfaunalAbundance and Master Species List** in query.
    Specify

-   desired Project,

-   desired Date,

-   AbundQualifier=Not Like \"L\*\",

-   Reason for Elimination=\"None\",

-   SCI_NM=Not Like \"Data lost\" AND Not Like \"Mollusca note\" AND Not
    Like \"No taxa found\".

> If not using all stations in the desired project and year, specify
> desired StationID list. If not using all field replicates, specify
> FieldReplicate=\"1\".

a)  Concatenate Project, Year, Station, and FieldReplicate, separated by
    > underscores, to create the Sample ID. Remove the space in \"Urban
    > Bays\".

b)  The fields required are Sample, SCI_NM, and Abundance.

```{=html}
<!-- -->
```
7.  **Check samples for usability:**

```{=html}
<!-- -->
```
a)  Match the Sample IDs in the benthos dataset to the sample IDs in the
    > Usability List (Sample, Usability(phylum), Usability(species/FG)).

b)  Keep only those samples for which Usability(species/FG)=TRUE.

```{=html}
<!-- -->
```
8.  **Combine with standardized taxa list** (SCI_NM, Standardized,
    TaxonR, Model):

```{=html}
<!-- -->
```
a)  Match the taxa names (SCI_NM) in the sample dataset to the SCI_NM
    taxa in the standardization list. It may be necessary to retroname
    taxa in the sample dataset to older names to match those in the
    standardization list.

b)  Keep the TaxonR names. (There may be more than one taxon with the
    same Standardized and TaxonR name.)

c)  Delete the taxa with the Standardized name \"delete\".

```{=html}
<!-- -->
```
9.  **Sum the abundances** for each sample in the sample dataset for
    each TaxonR name. This step must be performed in Minitab, Primer,
    or R. DO NOT USE A PIVOT TABLE IN EXCEL \-- a pivot table will
    provide incorrect sums.

10. **Prepare the summed, standardized data:**

```{=html}
<!-- -->
```
a)  Select the records of the summed abundances and the TaxonR name, for
    > just the 129 (or however many) taxa in the model.

b)  Pivot the selected, summed, standardized data to obtain a single
    > record for each sample containing the Sample ID and the summed
    > abundances for the 129 TaxonR taxa required for the model.

c)  Ensure that there are columns for all 129 taxa; add missing taxa, if
    > needed.

d)  Set all blank abundances to zero.

e)  Save as a csv file.

COMBINE ENVIRONMENTAL AND BENTHOS DATA

11. Match on the Sample IDs to concatenate the records to create a
    single record for each sample with the following fields: Sample,
    Depth, Penetration, Salinity, Temperature, Fines, Gravel, InorgC,
    and the abundances for the 129 TaxonR taxa required for the model.

12. Save as a csv file.
