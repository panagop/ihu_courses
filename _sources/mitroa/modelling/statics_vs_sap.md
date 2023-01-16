# Στατική vs SAP2000

## Διαφορές στα αποτελέσμτατα της επίλυσης με το χέρι και στην ανάλυση με το SAP2000

Αν επιχειρήσει να ελέγχει κάποιος τα αποτελέσματα μιας επίλυσης με το χέρι, με την κλασική στατική, εφαρμόζοντας μια από τις γνωστές μεθόδους (δυνάμεων, μετακινήσεων, cross, μητρώα κτλ), συγκρίνοντας με τα αποτελέσματα ενός προγράμματος πεπερασμένων στοιχείων, όπως το SAP2000, πιθανότατα θα παρατηρήσει κάποιες μικρές, όχι όμως ασήμαντες διαφορές.

Η διαφορές αυτές συχνά οφείλονται στην παραδοχή που γίνεται στην *κλασική* στατική να αγνοηθεί το έργο των τεμνουσών δυνάμεων. Έτσι, το μητρώο στιβαρότητας που χρησιμοποιείται για τα δομικά στοιχεία (όταν χρησιμοποιείται η μητρωική στατική) έχει την μορφή του παρακάτω σχήματος (ενώ οι ίδιες σχέσεις υπολογισμού εφαρμόζεται και στις άλλες μεθόδους).

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/k_matrix.png
---
width: 450px
align: center
name: k_matrix
---
```

Στην περίπτωση όμως που η επίλυση γίνει με ένα πρόγραμμα (π.χ. το SAP2000), η προεπιλεγμένη (default) επιλογή είναι να λαμβάνεται υπόψη και το έργο των τεμνουσών. Αυτό γίνεται χρησιμοποιώντας ένα διαφορετικό μητρώο στιβαρότητας (Timoshenko beam element [^timoshenko])

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/k_timoshenko.png
---
width: 600px
align: center
name: k_timoshenko
---
```

όπου $C=\dfrac{EI}{(1+φ)L^3}$ και $φ=\dfrac{12}{L^2} \dfrac{EI}{κGA}$

## Παράδειγμα με εναλλακτικές επιλύσεις

```{figure} https://raw.githubusercontent.com/panagop/ihu_courses/main/shared/images/modelling/shear_exclude/pi_foreas.png
---
width: 600px
align: center
name: pi_foreas
---
```

Ο φορέας του παρακάτω σχήματος θα επιλυθεί με εναλλακτικές προσεγγίσεις ώστε να γίνει εμφανής η επιρροή του συνυπολογισμού του έργου των τεμνουσών.

[^timoshenko]: Timoshenko, S.P.; Gere, J.M. Theory of Elastic Stability; McGraw-Hill: New York, USA, 1961
