# Functional traits

This module compiles and processes the functional trait data used in the analyses. The procedures implemented here follow the description provided in the Functional traits section of the Methods in the main manuscript and the Supplemental Methods section Functional traits of dung beetle species/morphospecies from Colombia.

Five functional traits were analyzed to describe key ecological strategies of dung beetle species: body size, front–rear leg ratio, nest guild, diet range, and diel activity. These traits represent morphological and behavioral axes through which species interact with dung resources, soil properties, and microclimatic conditions.

Forest–pasture conversion alters these environmental conditions by increasing habitat openness and thermal exposure, modifying soil structure and dung persistence, and simplifying vegetation cover. These changes act as environmental filters on traits associated with body size, locomotion and digging performance, dung relocation behavior, and activity period.

# Morphometric traits

Morphometric traits were derived from standardized measurements obtained from calibrated photographs processed using ImageJ. Nine measurements describing body dimensions and front and rear leg lengths were recorded for a subset of specimens (1–36 individuals per species).

From these measurements, two functional traits were calculated:

- Body size: Defined as body length × elytra width, used as a proxy for overall body mass and metabolic scaling.
- Front–rear leg ratio: Defined as (front femur length + front tibia length) / (rear femur length + rear tibia length + rear spur length). This ratio captures the relative investment in forelegs used for digging versus hind legs involved in locomotion and dung manipulation.

Individual measurements were averaged at the specimen level and subsequently aggregated to obtain species-level trait means, which were used in the functional diversity analyses.

The script Morphometrics_measurements.R compiles morphometric measurements obtained from multiple measurement campaigns and calculates species-level trait averages.

# Behavioral traits

Three additional functional traits were compiled from the literature and expert knowledge:

- Nest guild (tunneler, roller, or dweller)
- Diet range
- Diel activity (diurnal, nocturnal, or both)

Trait information was obtained through an extensive literature review and complemented with field-based expert knowledge. When species-specific information was unavailable, genus-level trait assignments were applied following established dung beetle trait syntheses.

# Trait dataset

The final trait dataset integrates morphometric measurements and behavioral traits at the species level and provides the functional trait matrix used for downstream analyses of functional diversity and environmental filtering.
