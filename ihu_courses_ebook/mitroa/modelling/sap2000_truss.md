# Επίλυση δικτυώματος με το SAP2000

![|600](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/mitroa/example9/ekfonisi.png)
![|600](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/mitroa/example9/ekfonisi.png")

## Γεωμετρία φορέα

Μια πιθανή επιλογή για τη γρήγορη εισαγωγή της γεωμετρίας του φορέα είναι ο συνδυασμός `2D Trusses` και στη συνέχεια `Slopped Truss`. Προφανώς θα μπορούσε να γίνει και άλλη επιλογή.

![|400](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_new_model.png)

![|500](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_2dtruss.png)

![|500](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_delete_horizontal.png)

Η οριζόντια ράβδος προφανώς πρέπει να αφαιρεθεί (επιλέγεται και στη συνέχεια `Delete`)

Χωρίς να είναι απαραίτητο, καλό είναι να δίνεται στο πρόγραμμα ότι η ανάλυση αφορά επίπεδο φορέα `Analyze` - `Set Analysis Options` - `Plane Frame`

![|500](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/sap_plane_frame.png)

## Υλικά

![|300](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_material.png)

### Συνθήκες στήριξης

Επιλέγονται οι κόμβοι στη βάση του φορέα και παγιώνονται οι μεταφορικοί βαθμοί ελευθερίας.

![|250](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_restraints.png)

## Διατομές

Η διατομή έχει εμβαδό $A=\dfrac{100}{\sqrt{2}}\text{mm}^2$ = $\dfrac{10^{-4}}{\sqrt{2}}\text{m}^2$ = $7.071\cdot 10^{-5}\text{m}^2$

Καθώς πρόκειται για δικτύωμα και μας ενδιαφέρουν μόνο οι αξονικοί βαθμοί ελευθερίας των ράβδων, οι ιδιότητες που σχετίζονται με την κάμψη, τη διάτμηση κτλ. μπορούν να αγνοηθούν (π.χ. ροπή αδράνειας).

Στο SAP2000 υπάρχει η δυνατότητα να ορίζεται *γενική διατομή* (`General Section`) στην οποία ο χρήστης θα εισάγει τις τιμές από τις ιδιότητες που τον ενδιαφέρουν. 

![|450](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_add_general.png)

![|450](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_general_section.png)

Οι διαστάσεις της γενικής διατομής, στην περίπτωσή μας δεν επηρεάζουν την επίλυση. Αντίστοιχα, στα `Section Properties` απαιτείται μόνο το εμβαδό της ράβδου.

![|450](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_property_data.png)

Ομοίως και στα modification factors, δίνεται 1 στην αξονική διεύθυνση, ενώ όλες οι υπόλοιπες ιδιότητες μηδενίζονται (το βήμα αυτό θα μπορούσε να παραληφθεί καθώς η σχετική πληροφορία ορίστηκε ήδη παραπάνω).

![|350](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_mod_factors.png)

### Ανάθεση των διατομών στις ράβδους του φορέα

![|500](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_frame_sections.png)

### Εισαγωγή των εσωτερικών αρθρώσεων στα άκρα των ράβδων

Επιλέγονται οι ράβδοι του δικτυώματος και με την επιλογή `Frame` - `Releases/Partial Fixity` ορίζεται ότι οι ροπές στα άκρα των ράβδων θα πρέπει να είναι μηδενικές.

![|450](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_releases.png)

![|400](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_releases2.png)

Όταν γίνει αυτό, η εποπτική απεικόνιση των εσωτερικών αρθρώσεων στο πρόγραμμα είναι σύμφωνα με το σχήμα που ακολουθεί.

![|500](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_releases_visual.png)

## Φορτία

### Load patterns/cases

![|450](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_load_patterns.png)

![|500](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_sap_load_cases.png)

### Εισαγωγή επικόμβιων φορτίων

![|500](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_assign_loads.png)

## Αποτελέσματα ανάλυσης

### Διάγραμμα αξονικών

![|600](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_results_axials.png)

### Μετακινήσεις κόμβου 2

![|600](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_results_joint_displ.png)

### Αντιδράσεις στήριξης

![|600](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_sap_results_reactions.png)