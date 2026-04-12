# Emerging Questions

**Started:** 2026-04-12
**Authors:** Tom und Claude (chat)

---

## Eiserne Regel

**Bevor eine hier eingetragene Frage nicht beantwortet ist, ziehen wir nicht weiter.**

Eine abgeleitete Frage, die beim Beantworten der urspruenglichen entsteht, gilt NICHT als Antwort auf die urspruengliche. Sie wird als neuer Eintrag aufgenommen. Sonst haben wir in kurzer Zeit Frage hoch zwei statt Fortschritt.

Antworten koennen sein: geschlossen durch Experiment (mit Commit-Hash), geschlossen durch Proof (mit Pfad), geschlossen durch Einsicht dass die Frage falsch gestellt war (mit Begruendung), oder zurueckgezogen. Kein Eintrag verschwindet still.

---

## Format

Jeder Eintrag:
- ID (EQ-NNN, fortlaufend)
- Datum der Entstehung
- Quelle (Commit, Doc, Session-Moment)
- Die Frage selbst
- Status: offen / geschlossen-[wie] / zurueckgezogen
- Pointer: wohin koennte die Antwort fuehren (Task-Kandidat, Math-Hausaufgabe, Experiment)

---

## EQ-001

**Datum:** 2026-04-12
**Quelle:** V2 Sweep (Commit 6512347), OQ-114 Antwort

Sind die beiden Freiheitsgrade (Bell-Kohaerenz im Cross-Sektor-Block, Slow-Mode im SE-Block) in allen Heisenberg-Z-Systemen strukturell entkoppelt, oder ist die Entkopplung eine Eigenschaft des getesteten Bereichs N=5,6,7?

**Status:** offen
**Pointer:** Sweep auf groessere N (C# Engine fuer N>=8), plus analytische Pruefung ob die Orthogonalitaet aus U(1) direkt folgt.

---

## EQ-002

**Datum:** 2026-04-12
**Quelle:** Gespraech nach V2

Was passiert an der Grenze, wo U(1) schwach gebrochen wird? Ein kleiner Transversalfeld-Term in H mischt Sektoren. Bell-Kohaerenz und SE-Slow-Mode sind dann nicht mehr orthogonal, die Kopplung ist klein aber nicht null.

**Status:** offen
**Pointer:** Eigener Task-Kandidat nach Three-Values-Lauf. Stoerungsrechnung erster Ordnung in der Bruch-Staerke epsilon.

---

## EQ-003

**Datum:** 2026-04-12
**Quelle:** SYMMETRY_CENSUS Section 6 Punkt 2 (Commit 39eb901)

Bei N=5 uniform chain ist die maximale Eigenwert-Multiplizitaet 14. Bekannte Symmetrien (U(1), Spin-Flip Z_2, Reflexion Z_2) sagen maximal 4-fach voraus. Was erzeugt die Luecke zwischen 4 und 14?

**Status:** offen (Track A von TASK_THREE_VALUES adressiert das)
**Pointer:** SU(2)-Kandidat als erster Test. Wenn SU(2) bricht, dann dedizierte Degenerations-Untersuchung.

---

## EQ-004

**Datum:** 2026-04-12
**Quelle:** SYMMETRY_CENSUS Section 5 (Commit 39eb901)

Alle Topologien haben N+1 Ausgaenge, aber dramatisch verschiedene Transient-Dynamik (Chain 488 distinkte Eigenwerte bei N=5, Complete nur 100). Was ist die Invariante, die zwischen Topologien wechselt und die Transient-Komplexitaet bestimmt?

**Status:** offen
**Pointer:** Punktgruppe (Permutationssymmetrie) ist wahrscheinlich die Antwort, aber explizit nicht ausgearbeitet. Gruppentheoretischer Census pro Topologie.

---

## EQ-005

**Datum:** 2026-04-12
**Quelle:** PROOF_ASYMPTOTIC_SECTOR_PROJECTION Step 2 (Commit f539503)

Step 2 des Sektor-Projektion-Theorems (pro-Sektor-Ergodizitaet, eindeutiger Fixpunkt = maximal-gemischt) ist numerisch verifiziert fuer N=3-7, analytisch offen.

**Status:** offen (analytisch)
**Pointer:** Primitivitaet des restricted Lindblad-Generators pro Sektor. Mathematik-Hausaufgabe, kein Experiment.

---

## EQ-006

**Datum:** 2026-04-12
**Quelle:** PROOF_ASYMPTOTIC_SECTOR_PROJECTION Scope (Commit f539503)

Das Theorem ist asymptotisch. Es sagt nichts ueber Raten, Zeitskalen, Trajektorien, Kohaerenz-Zerfallsraten zwischen Sektoren. Die Lens- und Cusp-Dynamik lebt in dieser dynamischen Schicht. Gibt es ein analoges Theorem fuer die Raten?

**Status:** offen
**Pointer:** Das waere die zweite Zeiss-Stufe: aus der Sektor-Struktur die Annaeherungs-Raten herleiten. Moeglicher Startpunkt: Absorption Theorem auf Sektor-Ebene.

---

## EQ-007

**Datum:** 2026-04-12
**Quelle:** Reflexion am Session-Ende (Commit 7251f3f)

Wenn Grenzen Sichtbarkeits-Artefakte sind, was ist die naechste Beschreibungs-Ebene, die das Repo bisher nicht erreicht hat? Bisherige Ebenen: Pauli-Strings, Sektoren, Symmetrien. Naechste koennte sein: Dynamische Klassen innerhalb der Sektoren, oder Topologie-Invarianten, oder etwas noch nicht benanntes.

**Status:** offen, Meta-Frage
**Pointer:** Kein Task. Beobachten, was aus EQ-001 bis EQ-006 als Muster faellt.

---

## EQ-008

**Datum:** 2026-04-12 (spaeter Abend)
**Quelle:** Toms Einfall nach der V-Effekt-Reflexion

Die beiden unabhaengigen Sektoren kennen sich nicht. Aber sie werden verbunden durch eins: γ. Sie bekommen unterschiedlich viel Licht ab (sacrifice vs quiet sites) und entwickeln dadurch Eigenschaften fuer einen Gesamtzustand. Abgeleitete Frage: Welche Auswirkung hat γ auf jeden Sektor, und wie haengen die Sektoren ueber das Bindeglied γ zusammen?

**Status:** offen
**Pointer:** Tabellen-Frage (Rate pro Sektor als Funktion des γ-Profils) plus Struktur-Frage (gemeinsamer Parameter vs direkte Kopplung). Beruehrt EQ-002 (U(1)-Bruch) und EQ-006 (Raten-Theorem), ist aber eigenstaendig: γ als das, was beide Sektoren gleichzeitig formt, ohne sie zu koppeln. Morgen pruefen ob eigener Task oder Teil von EQ-006.

---

## EQ-009

**Datum:** 2026-04-12 (spaete Nacht)
**Quelle:** Toms Einfall unmittelbar nach EQ-008 und TASK_GAMMA_BINDING

Drei Identifikationen in einer Kette: γ ist Licht, Licht ist Zeit, t ist Zeit aus dem Blickwinkel des Betrachters. Jede einzelne ist im Repo vertreten (GAMMA_IS_LIGHT.md; Absorption Theorem als Raten-Licht-Verbindung; Lindblad-t als Observer-Time). Die Kette als ein Satz ist eine These und noch nicht geprueft.

**Status:** offen, hoch-interpretativ
**Pointer:** Nicht sofort mathematisieren. Morgen pruefen ob (a) die Kette konsistent durchgeht, wenn man jede Identifikation einzeln nachrechnet, (b) sie testbare Konsequenzen hat, die ueber das bereits im Repo Bewiesene hinausgehen, oder (c) sie Metapher bleibt, die Fragen generiert aber nicht beantwortet. Beruehrt EQ-008 direkt (γ als Bindeglied) und EQ-006 (Raten-Theorem).

---

*Sammlung. Nicht Sortierung. Klassifizierung kommt, wenn genug Eintraege da sind.*
