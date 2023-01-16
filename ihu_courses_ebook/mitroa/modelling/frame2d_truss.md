# Επίλυση δικτυώματος με το Frame2d

![|600](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/mitroa/example9/ekfonisi.png)

To βασιμό μενού για την εισαγωγή των δεδομένων στο Frame2d είναι το `MODEL`

## Μονάδες

![|300](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_model_menu.png)

## Ιδιότητες δομικών στοιχείων

Επιλέγεται το σύστημα των μονάδων. Στο παράδειγμα θα χρησιμοποιηθεί το `Default Metric`

![|600](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_frame2d_units.png)

### Υλικό

Εισάγεται το νέο υλικό και παρατηρείται ότι η μόνη ιδιότητα που μας ενδιαφέρει είναι το μέτρο ελαστικότητας (δεν υπάρχει το μέτρο διάτμησης, ο λόγος του Poisson κτλ.).

![|450](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_frame2d_materials.png)

### Διατομές

Η διατομή έχει εμβαδό $A=\dfrac{100}{\sqrt{2}}\text{mm}^2$ = $\dfrac{1}{\sqrt{2}}\text{cm}^2$ = $0.7071\text{cm}^2$

Η ροπή αδράνειας δε μας χρειάζεται για το δικτύωμα, οπότε αφήνουμε τιμή $I=1.0\text{cm}^4$

![|450](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_frame2d_section.png)

## Γεωμετρία φορέα

### Κόμβοι

Οι κόμβοι δίνονται μέσω των συντεταγμένων τους και όχι με γραφικό τρόπο.
Οι συνθήκες στήριξης ορίζονται παγιώνοντας τους σχετικούς βαθμούς ελευθερίας. 

Στην περίπτωση που εξετάζεται υπάρχουν αρθρώσεις στη βάση του φορέα. Ωστόσο, σε πρώτη φάση επιλέγεται να παγιωθούν όλοι οι βαθμοί ελευθερίας (πακτώσεις) στους κόμβους 1(Α) και 3(C).

![|450](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_frame2d_nodes.png)

Στο σχήμα που ακολουθεί φαίνεται η τοποθέτηση και η ονομασία των κόμβων του φορέα.

![|450](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_frame2d_nodes2.png)

### Δομικά στοιχεία

Τα δομικά στοιχεία εισάγονται ενώνοντας τους κόμβους, επιλέγοντας το κατάλληλο υλικό και τη διατομή τους. 

Στην περίπτωση που υπάρχουν εσωτερικές αρθρώσεις στην αρχή ή στο τέλος ενός στοιχείου επιλέγεται το σχετικό άκρο στα `Hinges` (δεν υπάρχει άρθρωση στο φορέα που εξετάζεται).

Ωστόσο, επιλέγεται αυτό να μη γίνει ακόμα (θα γίνει στο επόμενο βήμα) και τα σχετικά πεδία αφήνονται κενά.

![|500](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_frame2d_elements.png)

Ως αποτέλεσμα αυτής της διαδικασίας, ο φορέας έχει προσομοιωθεί σύμφωνα με το σχήμα που ακολουθεί, με πακτώσεις στις στηρίξεις του, ενώ δεν έχουν ακόμα οριστεί οι αρθρώσεις στα άκρα των ράβδων του δικτυώματος.

![|450](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_frame2d_foreas_fixed.png)

## Φορτία

### Εισαγωγή επικόμβιων φορτίων

![|500](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_frame2d_nodal_loads.png)

## Μετατροπή του φορέα σε δικτύωμα

To Frame2D έχει τη δυνατότητα να μετατρέπει αυτόματα έναν πλαισιακό φορέα σε δικτύωμα με την επιλογή `TOOLS` - `Convert Model to Truss`

![|300](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/frame2d_convert_to_truss.png)

Το αποτέλεσμα αυτής της μετατροπής είναι να οριστούν αυτόματα από το πρόγραμμα οι αρθρώσεις στις στηρίξεις (αποδέσμευση της στροφής RZ) καθώς και η εσωτερική άρθρωση στο τέλος του στοιχείου 1(A) (δεν απαιτείται να μπει άλλη εσωτερική άρθρωση καθώς αυτό συνδέεται στον ίδιο κόμβο με το στοιχείο 2 ενώ υπάρχουν ήδη οι αρθρώσεις στις στηρίξεις).

![|450](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_frame2d_nodes_after_convert.png)

![|550](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_frame2d_elements_after_convert.png)

## Αποτελέσματα ανάλυσης

### Διάγραμμα αξονικών

![|550](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_frame2d_results_axials.png)

Αν τα διαγράμματα εμφανίζονται με περίεργο scaling (π.χ. είναι τεράστια), τότε με τα `+` και `-` κουμπιά, πάνω αριστερά, διορθώνεται το μέγεθος ώστε να εμφανίζονται λογικά.

### Μετακινήσεις κόμβου 2(B)

![|350](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_frame2d_results_joint_dipl.png)

### Αντιδράσεις στήριξης

![|550](https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/truss/truss_frame2d_results_reactions.png)

### Παρατηρήσεις

Τα αποτελέσματα ταυτίζονται με αυτά που προκύπτουν από το SAP2000 και είναι όμοια με αυτά που υπολογίζονται με την *κλασική στατική*.