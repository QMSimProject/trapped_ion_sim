#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    23.06.2014 10:36:54 CEST
# File:    gui.py

from qt_import import *
from .integration_routine import integrate_cf, integrate_cf_eff, plot

import numpy as np
import pickle

class Q2LaserWidget(QWidget):
    __laser_counter = 0
    def __init__(self, parent = None):
        self.__class__.__laser_counter += 1
        super(Q2LaserWidget, self).__init__(parent)
        self.init_ui()
    def close(self):
        self.__class__.__laser_counter -= 1
    def init_ui(self):
        #=================== widgets ===================
        self.name = QLineEdit(self)
        self.freq = QLineEdit(self)
        self.phase = QLineEdit(self)
        self.pattern = QLineEdit(self)
        #------------------- settings ------------------- 
        self.name.setPlaceholderText("name")
        self.freq.setPlaceholderText("freq")
        self.phase.setPlaceholderText("phase")
        self.pattern.setPlaceholderText("pattern")
        
        self.name.setStatusTip("name of the laser")
        self.freq.setStatusTip("frequency of the laser")
        self.phase.setStatusTip("phase of the laser")
        self.pattern.setStatusTip("firing pattern of the laser e.g: [[0, pi], [2*pi, 4.1]]")
        
        self.pattern.setMinimumWidth(130)
        
        self.name.setText("laser" + str(self.__class__.__laser_counter))
        #=================== layout ===================
        grid = QGridLayout()
        grid.addWidget(self.name   , 1, 1, 1, 1)
        grid.addWidget(self.freq   , 1, 2, 1, 1)
        grid.addWidget(self.phase  , 1, 3, 1, 1)
        grid.addWidget(self.pattern, 1, 4, 1, 1)
        
        self.setLayout(grid)
        
        self.show()

class Q2AtomWidget(QWidget):
    __atom_counter = 0
    def __init__(self, parent = None):
        self.__class__.__atom_counter += 1
        super(Q2AtomWidget, self).__init__(parent)
        self.init_ui()
    def close(self):
        self.__class__.__atom_counter -= 1
    
    def init_ui(self):
        #=================== widgets ===================
        self.name = QLineEdit(self)
        self.levels = QLineEdit(self)
        self.state = QLineEdit(self)
        
        #------------------- settings ------------------- 
        self.name.setPlaceholderText("name")
        self.levels.setPlaceholderText("energy levels")
        self.state.setPlaceholderText("start state")
        
        self.name.setStatusTip("name of the atom")
        self.levels.setStatusTip("energy levels e.g: [0, 10, 20]")
        self.state.setStatusTip("start state e.g: diag([1, 0, 0])")
        self.levels.setMinimumWidth(100)
        self.state.setMinimumWidth(100)
        
        self.name.setText("atom" + str(self.__class__.__atom_counter))
        #=================== layout ===================
        grid = QGridLayout()
        grid.addWidget(self.name   , 1, 1, 1, 1)
        grid.addWidget(self.levels , 2, 1, 1, 1)
        grid.addWidget(self.state  , 3, 1, 1, 1)
        
        self.setLayout(grid)
        self.show()
    
class Q2VibronWidget(QWidget):
    __vibron_counter = 0
    def __init__(self, parent = None):
        self.__class__.__vibron_counter += 1
        super(Q2VibronWidget, self).__init__(parent)
        self.init_ui()
    def close(self):
        self.__class__.__vibron_counter -= 1
    
    def init_ui(self):
        #=================== widgets ===================
        self.name = QLineEdit(self)
        self.freq = QLineEdit(self)
        self.state = QLineEdit(self)
        
        #------------------- settings ------------------- 
        self.name.setPlaceholderText("name")
        self.freq.setPlaceholderText("freq")
        self.state.setPlaceholderText("start state")
        
        self.name.setStatusTip("name of the vibron")
        self.freq.setStatusTip("frequency of the vibron")
        self.state.setStatusTip("start state e.g: diag([1, 0, 0])")
        self.state.setMinimumWidth(130)
        
        self.name.setText("vib" + str(self.__class__.__vibron_counter))
        #=================== layout ===================
        grid = QGridLayout()
        grid.addWidget(self.name   , 1, 1, 1, 1)
        grid.addWidget(self.freq , 1, 2, 1, 1)
        grid.addWidget(self.state  , 1, 3, 1, 1)
        
        self.setLayout(grid)
        
        self.show()

class Q2IntegrationWidget(QWidget):
    def __init__(self, parent = None):
        super(Q2IntegrationWidget, self).__init__(parent)
        
        #=================== widgets ===================
        self.interval = QLineEdit(self)
        self.measure = QLineEdit(self)
        self.cutoff = QLineEdit(self)
        
        #------------------- settings ------------------- 
        self.interval.setPlaceholderText("intervall")
        self.measure.setPlaceholderText("N")
        self.cutoff.setPlaceholderText("cutoff")
        
        self.interval.setStatusTip("integration intervall e.g: [0, 10*pi]")
        self.measure.setStatusTip("numbers of measurements")
        self.cutoff.setStatusTip("cutoff frequency for rwa")
        
        #=================== layout ===================
        grid = QGridLayout()
        
        grid.addWidget(self.interval     , 0, 0, 1, 1)
        grid.addWidget(self.measure      , 0, 1, 1, 1)
        grid.addWidget(self.cutoff       , 0, 2, 1, 1)
        self.setLayout(grid)
        
        self.show()

class Q2NumberWidget(QWidget):
    def __init__(self, parent = None):
        super(Q2NumberWidget, self).__init__(parent)
        
        #=================== widgets ===================
        self.atom_sb = QSpinBox(self)
        self.laser_sb = QSpinBox(self)
        self.vibron_sb = QSpinBox(self)
        
        self.atom_sb.setRange(1, 10)
        self.laser_sb.setRange(1, 20)
        self.vibron_sb.setRange(0, 10)
        #------------------- settings ------------------- 
        self.atom_sb.setStatusTip("number of atoms")
        self.laser_sb.setStatusTip("number of lasers")
        self.vibron_sb.setStatusTip("number of vibrons")
        
        #=================== layout ===================
        grid = QGridLayout()
        
        grid.addWidget(QLabel("Laser") , 1, 1, 1, 1)
        grid.addWidget(QLabel("Atom")  , 2, 1, 1, 1)
        grid.addWidget(QLabel("Vibron"), 3, 1, 1, 1)
        grid.addWidget(self.laser_sb  , 1, 2, 1, 1)
        grid.addWidget(self.atom_sb   , 2, 2, 1, 1)
        grid.addWidget(self.vibron_sb , 3, 2, 1, 1)
        self.setLayout(grid)
        
        self.show()

class Q2PlotWidget(QWidget):
    def __init__(self, parent = None):
        super(Q2PlotWidget, self).__init__(parent)
        grid = QGridLayout(self)
        
        grid.addWidget(QTextEdit(self), 1, 1, 1, 1)
        grid.addWidget(QLineEdit(self), 2, 1, 1, 1)
        
        self.setLayout(grid)
        self.show()
    
class Q2DisplayWidget(QMainWindow):
    def __init__(self, parent = None):
        
        super(Q2DisplayWidget, self).__init__(parent)
        self.laser_w = []
        self.atom_w = []
        self.vibron_w = []
        self.rabi_w = []
        self.eta_w = []
        self.init_ui()
        
    def init_ui(self):
        #=================== widgets ===================
        self.itg_w = Q2IntegrationWidget(self)
        self.num_w = Q2NumberWidget(self)
        
        #=================== connects ===================
        self.num_w.laser_sb.valueChanged.connect(self.changed_sb)
        self.num_w.atom_sb.valueChanged.connect(self.changed_sb)
        self.num_w.vibron_sb.valueChanged.connect(self.changed_sb)
        
        #=================== actions ===================
        exit_action = QAction(self.style().standardIcon(QStyle.SP_TitleBarCloseButton), "Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)

        load_action = QAction(self.style().standardIcon(QStyle.SP_DialogOpenButton), "Open", self)
        load_action.setShortcut("Ctrl+O")
        load_action.triggered.connect(self.load)
        
        save_action = QAction(self.style().standardIcon(QStyle.SP_DialogSaveButton), "Save as", self)
        save_action.setShortcut("Ctrl+S")
        save_action.triggered.connect(self.save)
        
        itg_action = QAction(self.style().standardIcon(QStyle.SP_MediaPlay), "Integrate", self)
        itg_action.setShortcut("Ctrl+P")
        itg_action.triggered.connect(self.integrate)
        
        #=================== menu bar ===================
        self.menuBar()
        menubar = self.menuBar()
        file_menu = menubar.addMenu("&File")
        file_menu.addAction(load_action)
        file_menu.addAction(save_action)
        file_menu.addAction(itg_action)
        file_menu.addAction(exit_action)
        
        #=================== tool bar ===================
        self.toolbar = self.addToolBar('File')
        #~ self.toolbar.addAction(exit_action)
        self.toolbar.addAction(load_action)
        self.toolbar.addAction(save_action)
        self.toolbar.addAction(itg_action)
        
        
        self.init_dynaminc_ui()
        
        self.statusBar()
        self.setWindowTitle("working title")
        #=================== layout ===================
    def init_dynaminc_ui(self):
        self.setup = QWidget(self)
        
        grid = QGridLayout()
        grid.addWidget(self.num_w     , 2, 1, 1, 1)
        
        #------------------- create lasers ------------------- 
        while self.num_w.laser_sb.value() < len(self.laser_w):
            x = self.laser_w.pop()
            x.close()
            del x
            y = self.rabi_w.pop()
            del y
        while self.num_w.laser_sb.value() > len(self.laser_w):
            self.laser_w.append(Q2LaserWidget(self))
            self.rabi_w.append([])
            for i in self.atom_w:
                self.rabi_w[-1].append(QLineEdit())
                self.rabi_w[-1][-1].setStatusTip("rabi frequency (pass a list e.g. [1, 2] if the rabi-freq are different for different transitions)")
        
        for l, l_i in zip(self.laser_w, range(len(self.laser_w))):
            grid.addWidget(l, 3 + l_i, 1, 1, 1)
        
        #------------------- create atoms ------------------- 
        while self.num_w.atom_sb.value() < len(self.atom_w):
            x = self.atom_w.pop()
            x.close()
            del x
            for r_row in self.rabi_w:
                y = r_row.pop()
                del y
            for e_row in self.eta_w:
                y = e_row.pop()
                del y
        while self.num_w.atom_sb.value() > len(self.atom_w):
            self.atom_w.append(Q2AtomWidget(self))
            for r_row in self.rabi_w:
                r_row.append(QLineEdit())
                r_row[-1].setStatusTip("rabi frequency (pass a list e.g. [1, 2] if the rabi-freq are different for different transitions)")
            for e_row in self.eta_w:
                e_row.append(QLineEdit())
                e_row[-1].setStatusTip("lamb-dicke parameter")
        
        for a, a_i in zip(self.atom_w, range(len(self.atom_w))):
            grid.addWidget(a, 2, 2 + a_i, 1, 1)
        
        #------------------- create vibrons ------------------- 
        while self.num_w.vibron_sb.value() < len(self.vibron_w):
            x = self.vibron_w.pop()
            x.close()
            del x
            y = self.eta_w.pop()
            del y
        while self.num_w.vibron_sb.value() > len(self.vibron_w):
            self.vibron_w.append(Q2VibronWidget(self))
            self.eta_w.append([])
            for i in self.atom_w:
                self.eta_w[-1].append(QLineEdit())
                self.eta_w[-1][-1].setStatusTip("lamb-dicke parameter")
        for v, v_i in zip(self.vibron_w, range(len(self.vibron_w))):
            grid.addWidget(v, 3 + v_i + len(self.laser_w), 1, 1, 1)
        
        #------------------- create rabi boxes ------------------- 
        for r_row, r_row_i in zip(self.rabi_w, range(len(self.rabi_w))):
            for r, r_i in zip(r_row, range(len(r_row))):
                grid.addWidget(r, 3 + r_row_i, 2 + r_i, 1, 1)
        #------------------- create eta boxes ------------------- 
        for e_row, e_row_i in zip(self.eta_w, range(len(self.eta_w))):
            for e, e_i in zip(e_row, range(len(e_row))):
                grid.addWidget(e, 3 + len(self.laser_w) + e_row_i, 2 + e_i, 1, 1)
        
        grid.addWidget(self.itg_w     , 3 + len(self.laser_w) + len(self.vibron_w), 1, 1, 1)
        self.setup.setLayout(grid)
        
        tab = QTabWidget(self)
        scroll = self.setup
        #~ scroll = QScrollArea(self)
        #~ scroll.setWidget(self.setup)
        tab.addTab(scroll, "Setup")
        tab.addTab(Q2PlotWidget(self), "Plot")
        self.setCentralWidget(tab)
        self.resize(0, 0)
        self.show()
        
    def changed_sb(self, val):
        self.init_dynaminc_ui()
    
    def parse(self):
        pi = np.pi
        diag = np.diag
        sqrt = np.sqrt
        
        cf = {}
        #=================== parse lasers ===================
        cf["laser"] = []
        for l in self.laser_w:
            name = l.name.text()
            freq = eval(l.freq.text())
            phase = eval(l.phase.text())
            pattern = eval(l.pattern.text())
            cf["laser"].append([name
                              , freq
                              , phase
                              , pattern
                               ])
        #=================== parse atoms ===================
        cf["atom"] = []
        for a in self.atom_w:
            name = a.name.text()
            levels = eval(a.levels.text())
            state = eval(a.state.text())
            assert len(state) == len(levels), 'dimension mismatch in ' + name
            cf["atom"].append([name
                             , levels
                             , state
                              ])
        #=================== parse vibrons ===================
        cf["vibron"] = []
        for v in self.vibron_w:
            name = v.name.text()
            freq = eval(v.freq.text())
            state = eval(v.state.text())
            cf["vibron"].append([name
                               , freq
                               , len(state)
                               , state
                               ])
        
        #=================== parse rabi ===================
        cf["rabi"] = []
        for r_row, r_row_i in zip(self.rabi_w, range(len(self.rabi_w))):
            temp = []
            for r, r_i in zip(r_row, range(len(r_row))):
                rabi = eval(r.text())
                n_lvl = len(cf["atom"][r_i][1])
                trans = n_lvl * (n_lvl - 1) / 2
                if isinstance(rabi, list):
                    assert len(rabi) == trans, 'dimension mismatch in rabi freq ' + str(rabi)
                else:
                    rabi = [rabi] * trans
                temp.append(rabi)
            
            cf["rabi"].append(temp)
        
        #=================== parse eta ===================
        cf["eta"] = []
        for i in self.atom_w:
            cf["eta"].append([])
        for e_row, e_row_i in zip(self.eta_w, range(len(self.eta_w))):
            for e, e_i in zip(e_row, range(len(e_row))):
                eta = eval(e.text())
                cf["eta"][e_i].append(eta)
        
        #=================== parse integration ===================
        cf["upper"] = eval(self.itg_w.interval.text())[1]
        cf["lower"] = eval(self.itg_w.interval.text())[0]
        cf["measure"] = eval(self.itg_w.measure.text())
        cf["cutoff"] = eval(self.itg_w.cutoff.text())
        
        return cf
    
    def save_parse(self):
        cf = {}
        #=================== parse lasers ===================
        cf["laser"] = []
        for l in self.laser_w:
            name = l.name.text()
            freq = l.freq.text()
            phase = l.phase.text()
            pattern = l.pattern.text()
            cf["laser"].append([name
                              , freq
                              , phase
                              , pattern
                               ])
        #=================== parse atoms ===================
        cf["atom"] = []
        for a in self.atom_w:
            name = a.name.text()
            levels = a.levels.text()
            state = a.state.text()
            cf["atom"].append([name
                             , levels
                             , state
                              ])
        #=================== parse vibrons ===================
        cf["vibron"] = []
        for v in self.vibron_w:
            name = v.name.text()
            freq = v.freq.text()
            state = v.state.text()
            cf["vibron"].append([name
                               , freq
                               , state
                               ])
        
        #=================== parse rabi ===================
        cf["rabi"] = []
        for r_row in self.rabi_w:
            temp = []
            for r in r_row:
                rabi = r.text()
                temp.append(rabi)
            
            cf["rabi"].append(temp)
        
        #=================== parse eta ===================
        cf["eta"] = []
        for i in self.atom_w:
            cf["eta"].append([])
        for e_row in self.eta_w:
            for e, e_i in zip(e_row, range(len(e_row))):
                eta = e.text()
                cf["eta"][e_i].append(eta)
        
        #=================== parse integration ===================
        cf["interval"] = self.itg_w.interval.text()
        cf["measure"] = self.itg_w.measure.text()
        cf["cutoff"] = self.itg_w.cutoff.text()
        
        return cf
        
    def load_parse(self, cf):
        self.num_w.laser_sb.setValue(len(cf["laser"]))
        self.num_w.atom_sb.setValue(len(cf["atom"]))
        self.num_w.vibron_sb.setValue(len(cf["vibron"]))
        #=================== parse lasers ===================
        for l_w, l in zip(self.laser_w, cf["laser"]):
            l_w.name.setText(l[0])
            l_w.freq.setText(l[1])
            l_w.phase.setText(l[2])
            l_w.pattern.setText(l[3])
        #=================== parse atoms ===================
        for a_w, a in zip(self.atom_w, cf["atom"]):
            a_w.name.setText(a[0])
            a_w.levels.setText(a[1])
            a_w.state.setText(a[2])
        #=================== parse vibrons ===================
        for v_w, v in zip(self.vibron_w, cf["vibron"]):
            v_w.name.setText(v[0])
            v_w.freq.setText(v[1])
            v_w.state.setText(v[2])
            
        #=================== parse rabi ===================
        for r_w_row, r_row in zip(self.rabi_w, cf["rabi"]):
            temp = []
            for r_w, r in zip(r_w_row, r_row):
                r_w.setText(r)
        
        #=================== parse eta ===================
        for e_w_row, e_w_row_i in zip(self.eta_w, range(len(self.eta_w))):
            for e_w, e_w_i in zip(e_w_row, range(len(e_w_row))):
                e_w.setText(cf["eta"][e_w_i][e_w_row_i])
        
        #=================== parse integration ===================
        self.itg_w.interval.setText(cf["interval"])
        self.itg_w.measure.setText(cf["measure"])
        self.itg_w.cutoff.setText(cf["cutoff"])
    
    def save(self):
        print("save")
        fileName = QFileDialog.getSaveFileName(self, "Open File", "~/", "Pickle Files (*.pickle)")
        file_name = fileName[0]
        if file_name == "":
            return
        if fileName[0].find(fileName[1].split("*")[-1][:-1]) == -1:
            file_name += fileName[1].split("*")[-1][:-1]
        
        data = self.save_parse()
        f = open(file_name, 'wb')
        pickle.dump(data, f)
        f.close()
        
        self.setWindowTitle(file_name.split("/")[-1])
    
    def load(self):
        print("load")
        fileName = QFileDialog.getOpenFileName(self, "Open File", "~/", "Pickle Files (*.pickle)")
        file_name = fileName[0]
        if file_name == "":
            return
        
        f = open(file_name, 'rb')
        data = pickle.load(f)
        f.close()
        self.load_parse(data)
        
        self.setWindowTitle(file_name.split("/")[-1])
    
    def integrate(self):
        
        canonical_form = self.parse()
        
        times, exp = integrate_cf_eff(canonical_form, method = "mesolve")
        plot(times, exp, canonical_form)
    
def run_gui():
    app = QApplication(sys.argv)
    w = Q2DisplayWidget()
    app.exec_()
